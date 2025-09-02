#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <pthread.h>
#include "common.h"
#include "qname_hash.h"

// Acces to bam / sam and bcf / vcf files.

// We want to hold connection data as an external pointer
// to the relevant R data structures:

struct bam_ptrs {
  samFile *sam;
  sam_hdr_t *header;
  hts_idx_t *index;
  // also include an iterator
  hts_itr_t *b_itr;
  // qname_hash is used to index query names to
  // alignment positions. It is a (string) ->
  // vector index; as such we need to make sure that
  // alignments to individual scaffolds are only indexed once.
  khash_t(sam_id_h) *qname_hash;
  khash_t(int) *indexed_targets;
};
// that will be held as an external pointer
// Which should be cleaned up

struct bcf_ptrs {
  vcfFile *vcf;  // can be bcf or vcf
  bcf_hdr_t *header;
  // the index only works with bcf; not sure about the iterator;
  hts_idx_t *index;
  hts_itr_t *b_itr;
};


static void finalise_bam_ptrs(SEXP ptr_r){
  //  Rprintf("finalising bam_ptrs\n");
  if(!R_ExternalPtrAddr(ptr_r))  // if already NULL, do nothing
    return;
  struct bam_ptrs *ptr = (struct bam_ptrs*)R_ExternalPtrAddr(ptr_r);
  // destroy the various indices
  if( ptr->header )
    sam_hdr_destroy( ptr->header );
  if( ptr->index )
    hts_idx_destroy( ptr->index );
  if( ptr->b_itr )
    hts_itr_destroy( ptr->b_itr );
  if(ptr->qname_hash)
    clear_sam_hash(ptr->qname_hash);
  if(ptr->indexed_targets)
    kh_destroy(int, ptr->indexed_targets);
  // It is _not_ safe to: 
  //  free( ptr->sam ); 
  // and there is no hts_file_destroy / sam_file_destroy
  sam_close( ptr->sam );
  free( ptr ); // ??
  R_ClearExternalPtr(ptr_r);
}

// It might be possible to finalise_bam_ptrs to handle bcf_ptrs, since the
// objects are quite similar. But it seems clearer to have a separate function
// call.
static void finalise_bcf_ptrs(SEXP ptr_r){
  if(!R_ExternalPtrAddr(ptr_r))  // if already NULL, do nothing
    return;
  struct bcf_ptrs *ptr = (struct bcf_ptrs*)R_ExternalPtrAddr(ptr_r);
  // destroy the various indices
  if( ptr->header )
    bcf_hdr_destroy( ptr->header );
  if( ptr->index )
    hts_idx_destroy( ptr->index );
  if( ptr->b_itr )
    hts_itr_destroy( ptr->b_itr );
  // It is _not_ safe to: 
  //  free( ptr->sam ); 
  // and there is no hts_file_destroy / sam_file_destroy
  free( ptr ); // ??
  R_ClearExternalPtr(ptr_r);
}


// Convenience function that makes a char sequence out of
// the sequence in bam1_t.data
// It would be better for this to take a reference to a buffer
// that can be reallocated if needed. As it is we keep calling
// malloc and free which is rather inefficient.
char *bam_seq(bam1_t *bam){
  uint8_t *seq_data = bam_get_seq(bam);
  // the sequence may not be stored in the bam file. In this case
  // it is not clear as to what bam_get_seq will return, but a value of
  // 0 is would not be strange.
  int32_t seq_l = bam->core.l_qseq;
  if(!seq_l)
    return(0);
  char *seq = malloc( seq_l + 1 );
  seq[seq_l] = 0;
  for(int i=0; i < seq_l; ++i)
    seq[i] = nuc_encoding[ bam_seqi(seq_data, i) ];
  return(seq);
}

// Avoid malloc if possible:
// but note that whatever we do, we end up copying data too much
// Note that we could probably use, a single call of `nibble2base`
// to get the sequence. But it doesn't look like nibble2base is exported
size_t bam_seq_p(bam1_t *bam, char **seq, size_t seq_l){
  if(seq_l < 1 + bam->core.l_qseq){
    seq_l = 1 + bam->core.l_qseq;
    *seq = realloc( *seq, 1 + seq_l );
  }
  uint8_t *seq_data = bam_get_seq(bam);
  (*seq)[bam->core.l_qseq] = 0;
  for(int i=0; i < bam->core.l_qseq; ++i)
    (*seq)[i] = nuc_encoding[ bam_seqi(seq_data, i) ];
  return(seq_l);
}

size_t bam_qual_p(bam1_t *bam, char **qual, size_t qual_l){
  if(qual_l < 1 + bam->core.l_qseq){
    qual_l = 1 + bam->core.l_qseq;
    *qual = realloc( *qual, qual_l );
  }
  (*qual)[ bam->core.l_qseq ] = 0;
  memcpy( *qual, bam_get_qual(bam), bam->core.l_qseq );
  return(qual_l);
}

size_t bam_aux_p(bam1_t *bam, char **aux, size_t buf_l){
  size_t aux_l = bam_get_l_aux(bam);
  if(buf_l < 1 + aux_l){
    buf_l = 1 + aux_l;
    *aux = realloc( *aux, buf_l );
  }
  (*aux)[ aux_l ] = 0;
  memcpy( *aux, bam_get_aux(bam), aux_l );
  return(buf_l);
}

// Appends a character representation of an auxiliary field to a 
// kstring_t object. kstring_t is an htslib struct.
// returns: the length of the resulting 0 terminated string.
// Note that use of kstring_t means a _lot_ of calls to `realloc`
// and that we have to clear the kstring_t memory after use.
// The code here is taken more or less directly from 
// the relevant section of sam_format1_append()
size_t bam_aux_string(bam1_t *b, kstring_t *str){
  int r = 0;
  uint8_t *s, *end;
  //  const bam1_core_t *c = &b->core;

  s = bam_get_aux(b); // aux
  end = b->data + b->l_data;
  
  while (end - s >= 4) {
    // only add a tab character if the string has length
    if(str->l)
      r |= kputc_('\t', str);
    if ((s = (uint8_t *)sam_format_aux1(s, s[2], s+3, end, str)) == NULL)
      goto bad_aux;
  }
  r |= kputsn("", 0, str); // nul terminate
  // here we ignore that, but it is something that should be handled.
  return(str->l);
 bad_aux:
  /* warning("Corrupted aux data for read %.*s", */
  /* 	  b->core.l_qname, bam_get_qname(b)); */
  //  errno = EINVAL;
  return -1;
}

// helper function for parse_MM_string
// determines the modification code; sets the modification indicated
// as an unsigned integer made by combining up to four modifications.
// returns a pointer to where the data starts. That is, after
// the first comma. Also sets the char skip_info to "?" or "."
// of ',' if not defined.
//
// *data should point to the first character of the modification code
//  that comes after the +/- indication.
// returns 0 on 
uint8_t *get_ml_code(uint8_t *data, uint32_t *mod_code, int *mod_n, char *skip_info){
  *mod_code = 0;
  *mod_n = 0;
  *skip_info = ',';
  size_t i=0;
  while(data[i] && data[i] != '?' && data[i] != '.' && data[i] != ','){
    ++i;
  }
  if(data[i] == 0 || i == 0){
    // this will be true at the end of the MM entry as it should end with ;
    // hence do not warn here:
    // warning("Found end of string in MM header");
    return(0);
  }
  if(i > 4){
    //    warning("No more than four modifications per position allowed");
    return(0);
  }
  *skip_info = data[i];
  if(data[0] >= '0' && data[0] <= '9'){
    *mod_n = 1;
    *mod_code = atoi( (const char*)data );
  }else{
    *mod_n = i;
    for(size_t j=0; j < i; ++j)
      *mod_code |= (((uint32_t)data[j]) << (j * 8));
  }
  if(data[i] == ',')
    return( data + i + 1 );
  return( data + i + 2 );
}

// translate query to reference positions. Requires an i_matrix of
// query operations; 
// Returns transformed query position (1 based) on alignment, -1 on no alignment
// to the one used that contains the query position.
// q_pos: the query position. This is assumed to be 1-based because of
//        the behaviour of parse_MM_string
// r_pos: A pointer whose value should be set to the reference position
// ops:   A pointer to an i_matrix table containing cigar positions
//        should also use 1 based offsets.
// beg:   The first column of the current alignment in the ops table
// end:   The index of the last column of the current alignment + 1
// col_i: The current column being parsed. This should be incremented
//        or decremented depending on the circumstance.
// is_fwd: 0 if the read has been reverse complemented. In that case
//         we have to go backwards.
// qlen:  The length of the query sequence. Needed for reverse mapping.
// in_clipped: what to do if the query position is in a clipped region
//             0: return -1
//             1: return the inferred position
//             2: return the position of the adjacent match operation.
int query_to_ref(int32_t q_pos, int *r_pos, struct i_matrix *ops,
		 size_t beg, size_t end, size_t *col_i, int is_fwd, int32_t q_len,
		 int in_clipped){
  // Positions obtained from from parse_MM_string are 1 based due
  // to how they are defined; the positions in the ops table used to be
  // 0-based, but are now 1 based to fit in with R semantics.
  // default mapping:
  *r_pos = -1;
  if(beg == end || *col_i == end || *col_i < beg)
    return(-1);
  // the row offsets of the ops table; We might want to check
  // the structure, but this should be considered a priviledged
  // call.
  int op = 1;
  int r0 = 3; 
  int q0 = 4;
  int r1 = 5;
  int q1 = 6;
  if(is_fwd){
    while( *col_i < end && *col_i < ops->col){
      int *coords = ops->data + (*col_i) * ops->nrow;
      if( coords[q0] <= q_pos && coords[q1] > q_pos ){
	// if the op type is M, (i.e. 0) we have a match.
	// else consider we have clipped region:
	if( coords[op] == 0 ){
	  *r_pos = coords[r0] + (q_pos - coords[q0]);
	  return(q_pos);
	}
	switch(in_clipped){
	case 1:
	  *r_pos = coords[r0] + (q_pos - coords[q0]);
	  break;
	case 2:
	  *r_pos = (*col_i) > beg ? coords[r0] : coords[r1];
	  break;
	default:
	  *r_pos = -1;
	}
	return(q_pos); // q_pos is within a cigar op
      }
      (*col_i)++;
    }
    return(-1); // q_pos is outside of the cigar ops: this shouldn't happen
  }
  // otherwise we are mapping in the reverse direction.
  // We note that we might be able merge the two parts here, as the only difference is
  // the loop condition and the direction of change. But lets try separately first.
  q_pos = 1 + q_len - q_pos;
  if(q_pos < 0){
    //    warning("Obtained a negative query position: %d for q_len %d", q_pos, q_len);
    return(-1);
  }
  while( *col_i >= beg ){
    if(*col_i >= ops->col){
      //      warning("Inappropriate value of *col_i (negative?): %ld", *col_i);
      return(-1);
    }
    int *coords = ops->data + (*col_i) * ops->nrow;
    if( coords[q0] <= q_pos && coords[q1] > q_pos ){
      // if the op type is M, (i.e. 0) we have a match.
      // otherwise we don't have a match.
      if( coords[op] == 0 ){
	// do not add 1 here as it's already added to q_pos above.
	*r_pos = coords[r0] + (q_pos - coords[q0]);
	return(q_pos);
      }
      switch(in_clipped){
      case 1:
	*r_pos = coords[r0] + (q_pos - coords[q0]);
	break;
      case 2:
	*r_pos = (*col_i) > beg ? coords[r0] : coords[r1];
	break;
      default:
	*r_pos = -1;
      }
      return(q_pos);
    }
    (*col_i)--;
  }
  return(-1);
}

// 
// This should be replaced with a function allowing the extraction
// of arbitrary auxiliary data; but I do not have sufficient time
// right now.
// AS: alignment score
// NM: edit distance
// NH: number of reported alignments
// IH: query hit total counts
// AS and NM appear common; NH and IH seem useful. This is experimental
// and should be replaced at some point later on.. 
void extract_int_aux_values(bam1_t *b, int *AS, int *NM, int *NH, int *IH){
  // actually I should use
  // bam_aux_get() and bam_aux2i() here
  const char as[4][2] = {"AS", "NM", "NH", "IH"};
  int *ret[4] = {AS, NM, NH, IH};
  for(int i=0; i < 4; ++i){
    *ret[i] = R_NaInt;
    uint8_t *s = bam_aux_get(b, as[i]);
    if(!s)
      continue;
    // NOTE: bam_aux2i returns a int64_t; but R does not have int64
    //       and it is unlikely that we will need it for the tags used
    //       here.
    // bam_aux2i returns 0 if there was an error; detecting the error
    // requires checking the errno. 
    errno = 0;
    int64_t v = bam_aux2i(s);
    *ret[i] = (errno != EINVAL) ? (int)v : R_NaInt;
  }
}

// Returns expected number of erros in query sequence based on
// a Poisson approximation of the Poisson binomial distribution
// The return value is expressed as (256 * expected_number) in order
// to provide non-integral numbers. 256 was chosen because of its
// similarity with ML (modification likelihood) encoding and because
// of the low likelihood of overflow (2^31 / 2^8 = 2^23; i.e. an
// overflow will not happen in R unless there are 2^23 or more expected
// errors; that's 8 million errors.)
//
// The length argument uses int32_t becaus this is the type held by a
// bam1_t object
int32_t estimate_read_errors(const unsigned char* qual, int32_t ql){
  const double mult = 256;
  double err_n = 0;
  for(int32_t i=0; i < ql; ++i){
    err_n += pow(10, (-((double)qual[i]) / 10.0));
  }
  return((int32_t)( err_n * mult ));
}

// parse MM AND ML strings; fill a  matrix with values;
// 
// This function assumes that all information will be contained within a single
// MM and a single ML AUX field. I do not know if this is what the standard
// specifies. This should be checked.
//
// s is 0 or a position in the bam_aux_string (ideally where MM starts)
// if 0, it will be set to the beginning of the AU data;
// I should probably remove this argument as it is possible for ML to precede MM
// in which case specifying s is problematic.
//
// To define positions in the reference we need to parse the cigar; we could do that
// directly here, but the complete cigar needs to be parsed first for reverse sequences
// That means that we might as well use the existing code that does this, and take
// a table of cigar operations with query and reference positions already defined.
// cig_ops:   a pointer to an i_matrix table filled by cigar_to_table()
// ops_beg and ops_end: the range of rows giving the cigar operations
// for the alignment specified by b. Note that ops_end is the index of
// the first operation for the next query; hence ops_end - ops_beg gives
// the numer of operations; note that this can be 0. This can happen for unaligned
// query sequences and must be checked.
// query_to_ref(cigar, cigar_pos)
// that moves to the next cigar position when needed to might make the most sense

// NOTES: I should be able to remove the cig_ops, ops_begin, ops_end arguments
// as I can use b->core.n_cigar instead. However, this function calls query_to_ref
// and I would need to refactor that function as well.
// ref_seq can b 0.
void parse_MM_string(bam1_t *b, uint8_t *s, struct i_matrix *mod_data, int al_i,
		     struct i_matrix *cig_ops, size_t ops_beg, size_t ops_end,
		     const char *ref_seq, size_t ref_seq_l){
  // mod_data should have 6 rows; these are
  // al_i, query_pos, mod_code, mod_n, mod_likelihood, ref_pos
  // This is very wasteful; we could get rid of al_i and return as elements of a list instead
  // but that would just move the problem downstream.
  //
  if(mod_data->nrow != MM_INFO_RN){
    return;
  }
  if(ops_beg >= ops_end || ops_end > cig_ops->col){
    return;
  }
  int mod_data_begin = mod_data->col; 
  // return if hard clipped. With hard clipping we have
  // no means of identifying the locations in the sequence.
  uint32_t *cigar = bam_get_cigar(b);
  // Remove any with hard clipping set.
  if(bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP || cig_ops->data[ (ops_end-1) * cig_ops->nrow + 1 ] == BAM_CHARD_CLIP){
    return;
  }
  int is_fwd = (b->core.flag & 16) ? 0 : 1;
  size_t current_op = is_fwd ? ops_beg : ops_end - 1;
  // MM tags take the following pattern: "MM:Z:C+h?,
  // with the ':' omitted in the binary data.
  // bam_aux_first, bam_aux_next give pointers to the beginning of the auxiliary
  // string.
  // the key (MM), type (Z) and tag can be obtained from: s - 2, *s and s+1
  // that is to say, that the complete field starts at s -2
  // since we need the key and type we simply subtract 2 from s.
  // but we cannot subtract before checking for null
  if(s == 0)
    s = bam_aux_first(b);
  uint8_t *end = b->data + b->l_data;
  short MM = *((short*)"MM");
  short ML = *((short*)"ML");
  // mm and ml hold the positions of the MM and ML auxiliary data.
  uint8_t *mm = 0;
  uint8_t *ml = 0;
  while(s){ 
    if(*(short*)(s-2) == MM)
      mm = s-2;
    if(*(short*)(s-2) == ML)
      ml = s-2;
    s = bam_aux_next(b, s);
  }
  if(mm == 0){
    // either some error, or no data; note that ml, is optional.
    // we should consider checking errno; if EINVAL, then data is corrupt.
    return;
  }
  // The modification code should start at mm+5; this must not be past the end of the
  // auxiliary data since we directly take the address.
  if(mm + 5 >= end){
    return;
  }
  // the type for MM should be Z. We don't know how many positions we are going
  // to parse; 
  // the type for ML should be B, with C, and a 4 bit count 
  // The mm : "MMZ" Followed by:
  // ([ACGTUN][-+]([a-z]+|[0-9]+)[.?]?(,[0-9]+)*;)*
  if(mm[2] != 'Z'){
    return;
  }
  char nuc = mm[3];
  //  char strand = mm[4];
  // determine how many modifications are related to each entry
  uint32_t mod_code;
  char skip_info;
  int mod_n;
  uint8_t *data = get_ml_code(mm + 5, &mod_code, &mod_n, &skip_info);
  uint8_t *seq_data = bam_get_seq(b);
  uint8_t *seq_qual = bam_get_qual(b);
  if(seq_data == 0 || seq_qual == 0){
    return;
  }
  uint32_t nuc_info = 0; // will hold, ref, query base and query quality
  int32_t seq_i = 0;
  int32_t qpos = 0;
  int32_t seql = b->core.l_qseq;
  uint32_t mod_pos_count = 0;
  // in theory it's possible for the pointer to become NULL, but
  // that would have to be after exceeding 2^64 - 1. Which would
  // mean that it's invalid anyway as we won't have that much memory.
  while(data && data < end && *data >= '0' && *data <= '9'){
    char base = 0;
    int nb = atoi((const char*)data);
    mod_pos_count++;
    // then process the query sequence to find the query position
    int base_count = 0;
    int r_pos = -1;
    if(nuc == 'N'){
      seq_i = seq_i + nb;
    }else{
      while(seq_i < seql && base_count <= nb){
	base = is_fwd ? nuc_encoding[ bam_seqi(seq_data, seq_i) ] : nuc_encoding_rc[ bam_seqi(seq_data, seql-(seq_i+1)) ];
	if( base == nuc )
	  ++base_count;
	++seq_i;
      }
    }
    // seq_i is one based; this means that it should be OK for it to be equal
    // to seq_l
    if(seq_i > seql){
      break;
    }
    // At this point the base at seq_i should be nuc; and should
    // be the one for which we have modification data. 
    // add a column to a table with the following rows:
    // alignment index, query position (seq_i), modification code, mod_n (number of modifications, up to 4), mod_likelihood, ref_pos
    // obtain the refeference position.
    int in_clipped = 0; // set to -1 if operation is not in an M operation
    qpos = query_to_ref( seq_i, &r_pos, cig_ops, ops_beg, ops_end, &current_op, is_fwd, seql, in_clipped );
    nuc_info = 0;
    char fwd_base = nuc_encoding[ bam_seqi(seq_data, qpos-1) ];
    if(qpos > 0){
      nuc_info |= (((uint32_t)base) << 24);
      if(seq_qual && (size_t)qpos <= seql)
	nuc_info |= (((uint32_t)seq_qual[ qpos - 1 ]) << 16);
      if(ref_seq && (size_t)r_pos <= ref_seq_l){
	nuc_info |= (((uint32_t)ref_seq[ r_pos - 1 ]) << 8);
	nuc_info |= (ref_seq[r_pos-1] == fwd_base);
      }
      // if the base is fwd, set the second bit:
      nuc_info |= (is_fwd << 1);
    }
    push_column( mod_data, (int[]){al_i, (int)qpos, (int)mod_code, mod_n, 0, r_pos, nuc_info} );
    // find the comma or the semicolon
    while(*data && *data >= '0' && *data <= '9' )
      ++data;
    if(*data == ';'){
      if(data[1] == 0)
	break;
      nuc = data[1];
      //      strand = data[2];
      data = get_ml_code(data+3, &mod_code, &mod_n, &skip_info);
      current_op = is_fwd ? ops_beg : ops_end - 1;
      // here data can be 0 if the ';' indicates the end of the MM field.
      seq_i = 0;
      continue;
    }
    ++data;
  }
  // if ml is 0, then we can't do much more
  if(ml == 0)
    return;
  // mm is encoded as ML:B:C, that is as single bytes; we may have more than one if
  // mod_n > 1. Should not be more than 4.
  if(ml[2] != 'B'){
    return;
  }
  if(ml[3] != 'C'){
    return;
  }
  // An auxiliary with a B type has the following structure:
  // [tag:2 bytes][B: 1 byte][ c|C|s|S|i|I|f 1 byte ][ count: 4 bytes (int) ]
  // Hence we get the number of entries at ml + 4
  uint32_t ml_count = *((int32_t*)(ml + 4));
  //  if(ml_count != mod_pos_count)
    /* warning("parse_MM: ml_count (%d) does not equal mod_pos_count (%d) for %s\n(end - ml = %ld)", */
    /* 	    ml_count, mod_pos_count, bam_get_qname(b), end-ml); */
  // And the quality (likelihood of modification data) starts at ml + 8
  uint8_t *qual = ml + 8;
  uint32_t i=0;
  int col = mod_data_begin;
  while(i < ml_count){
    if(col >= mod_data->col){
      //      warning("exceeded the number of columns in mod_data");
      break;
    }
    int *col_data = mod_data->data + col * mod_data->nrow;
    for(int j=0; j < col_data[3]; ++j){
      col_data[ 4 ] |= (((uint32_t)(qual[i])) << (j*8));
      ++i;
    }
    ++col;
  }
  // At this point we have certain expectations. For example
  // col should be equal to mod_data->col.
  /* if( col != mod_data->col ){ */
  /*   warning("col %d is not equal to mod_data->col %ld for %s", col, mod_data->col, bam_get_qname(b)); */
  /* } */
  return;
}


// Parse the cigar data and extend:
// A cigar operations table containing the operations, their lengths and
//     the query and reference positions (i.e. how they align)
// Sets q_length the query length calculated from the cigar string; this can differ
//     from that given in the bam1_t data struct due to hard clipping.
//     Similarly, the beginning and end positions of the original query can differ
//     So the function also sets q_beg and q_end. These refer to the query
//     prior to any hard clipping.
//
// Set the value of ops_beg, and ops_end giving the first and end columns of
// al_coord for the given alignment. These can then be used when parsing
// MM and ML tags in order to get reference positions.
// 
// ops_end is defined as 1 + the last index; this is so that ops_beg == ops_end when
// no operations were parsed. This can happen for unaligned reads.
// 
// Note that this implies that parsing MM / ML requires that the cigar is
// parsed.
//
void cigar_to_table(bam1_t *al, int al_i, int *r_pos,
		    struct i_matrix *al_coord,
		    int *qcig_length, int *q_beg, int *q_end,
		    size_t *ops_beg, size_t *ops_end,
		    struct cigar_parse_options *opts, int do_push){
  // Only add to al_coord / diff / mm_info tables if do_push is TRUE
  
  // use a 1-based coordinate system instead of a 0-based one
  // for consistency with R and aligned_regions
  *r_pos = 1 + (int)al->core.pos;
  // q_pos refers to the actual position in the original query string (in
  // 1-based coordinates) and will be returned as part of the ops table
  // If the sequence has been hard clipped, then we will not have a record
  // of it in this alignment record. In this case we need to adjust for the
  // hard clipping to get the corect position.
  int q_pos = 1;
  int h_clip = 0;
  *qcig_length = 0;
  *q_beg = 1;
  *q_end = 0;

  uint32_t *cigar = bam_get_cigar(al);
  uint32_t cigar_l = al->core.n_cigar;
  if(cigar_l < 1)
    return;
  // if we have clipping set q_beg appropriately:
  // 5 XOR [4, 5] --> 0, 1; this tests for BAM_CSOFT_CLIP and BAM_CHARD_CLIP
  *q_beg += (bam_cigar_op(cigar[0]) ^ 5) < 2 ? bam_cigar_oplen( cigar[0] ) : 0;
  if(bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP){
    // If the first cigar op is hard clip, then it _should_ guarantee the presence
    // of a second cigar op; hence we should not need to check the cigar length here.
    *q_beg += bam_cigar_op(cigar[1]) == BAM_CSOFT_CLIP ? bam_cigar_oplen(cigar[1]) : 0;
    h_clip = bam_cigar_oplen( cigar[0] );
  }
  
  *ops_beg = al_coord->col;
  int r0, q0;
  for(uint32_t i=0; i < cigar_l; ++i){
    // if operation is H also increment the qcig_length
    uint32_t type = bam_cigar_type(cigar[i]);
    uint32_t op = bam_cigar_op(cigar[i]);
    uint32_t op_length = bam_cigar_oplen( cigar[i] );
    r0 = *r_pos;
    q0 = q_pos;
    
    if( type & 1 || op == BAM_CHARD_CLIP ){
      *qcig_length += op_length;
      q_pos += op_length;
    }
    //q_pos += (type & 1) || op == BAM_CHARD_CLIP ? op_length : 0;
    *r_pos += type & 2 ? op_length : 0;
    if(op == BAM_CMATCH)
      *q_end = q_pos;
    if(do_push)
      push_column( al_coord, (int[]){ al_i, op, type, r0, q0, *r_pos, q_pos, op_length } );
    if((opts->depth || opts->diff) && op == BAM_CMATCH){
      // the r0 and q0 coordinates are 1 based; we need to adjust for this
      for(uint32_t j=0; j < op_length; ++j){
	size_t r = (r0 + j) - 1;
	size_t q = (q0 + j) - (1 + h_clip);
	size_t d_o = r - opts->depth_begin;
	// negative coordinates should overflow and hence be larger than the buffers
	// unless these are extremely large (close to the limit of size_t)
	if(opts->depth && d_o < opts->depth_l) ++opts->depth[d_o];
	if(!do_push)
	  continue;
	// Here I assume that if qseq_l is set then qseq will also be set
	// I leave qqual as an option as it will just mean quality = 0
	if(opts->diff && opts->qseq && q < opts->qseq_l
	   && r < opts->rseq_l && (opts->qseq[q] != opts->rseq[r])){
	  uint32_t qq = opts->qqual ? (opts->qqual[q] << 16) : 0;
	  int nuc_info = (int)( ((uint32_t)opts->rseq[r] << 8)) | ((uint32_t)opts->qseq[q]) | qq;
	  push_column( opts->diff, (int[]){ al_i, r+1, q+1, nuc_info });
	}
      }
    }
    if(opts->i_depth && op == BAM_CREF_SKIP && (op_length <= opts->max_intron_length)){
      for(int j=0; j < op_length; ++j){
	size_t d_o = r0 + j - (1 + opts->depth_begin);
	if(d_o < opts->depth_l)
	  ++opts->i_depth[d_o];
      }
    }
  }
  if(!do_push)
    return;
  *ops_end = al_coord->col;
  if(opts->mm_info){
    parse_MM_string(al, 0, opts->mm_info, al_i, al_coord, *ops_beg, *ops_end, opts->rseq, opts->rseq_l);
  }
}

// A convenience function; creates a VECSXP containing a STRSXP in the first
// element that can be used with setAttrib for R_DimNamesSymobl when the number 
// of rows is known.
SEXP mk_rownames(const char **names, size_t n){
  SEXP dims = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT( dims, 0, allocVector(STRSXP, n));
  for(size_t i=0; i < n; ++i)
    SET_STRING_ELT( VECTOR_ELT(dims, 0), i, mkChar(names[i]));
  UNPROTECT(1);
  return(dims);
}

SEXP mk_strsxp(const char **words, size_t n){
  SEXP words_r = PROTECT(allocVector(STRSXP, n));
  for(size_t i=0; i < n; ++i)
    SET_STRING_ELT(words_r, i, mkChar(words[i]));
  UNPROTECT(1);
  return(words_r);
}

// if index is not found, then warn and set the index to 0
// but still expect a character of some sort..
SEXP load_bam(SEXP bam_file_r, SEXP index_file_r){
  if(TYPEOF(bam_file_r) != STRSXP || length(bam_file_r) != 1)
    error("bam_file_r should contain a single string giving the bam file name");
  if(TYPEOF(index_file_r) != STRSXP || length(index_file_r) != 1)
    error("index_file_r should contain a single string giving the name of the bam index file");
  const char *bam_file = CHAR(STRING_ELT( bam_file_r, 0 ));
  const char *index_file = CHAR(STRING_ELT( index_file_r, 0 ));
  // make a bam_ptrs and load files and indices before returning these..
  samFile *sam = sam_open( bam_file, "r" );
  if(!sam)
    error("Unable to open specified bam file");
  hts_idx_t *index = sam_index_load(sam, index_file);
  if(!index){
    warning("No index file set; operations requiring index will not be available");
      //    error("Unable to open index");
  }
  sam_hdr_t *header = sam_hdr_read(sam);
  // Do I need to allocate this on the stack?
  struct bam_ptrs *ptr = malloc( sizeof(struct bam_ptrs) );
  ptr->sam = sam;
  ptr->header = header;
  ptr->index = index;
  ptr->b_itr = 0;  // to set an iterator we have to make a range query..
  ptr->qname_hash = 0; // holds a hash of query name -> list of alignment positions.
  ptr->indexed_targets = 0;
  // create an external pointer with suitable tag and protect?
  SEXP tag = PROTECT( allocVector(STRSXP, 1) );
  SET_STRING_ELT( tag, 0, mkChar("bam_ptrs") );
  // Would suggest to make prote R_NilValue
  // As I don't understand the purpose of the prot field
  SEXP prot = PROTECT( allocVector(STRSXP, 2) );
  SET_STRING_ELT( prot, 0, STRING_ELT(bam_file_r, 0));
  SET_STRING_ELT( prot, 1, STRING_ELT(index_file_r, 0));
  SEXP ptr_r = PROTECT( R_MakeExternalPtr(ptr, tag, prot) );
  R_RegisterCFinalizerEx(ptr_r, finalise_bam_ptrs, TRUE);
  UNPROTECT(3);
  return(ptr_r);
}

// This is almost identical to load_bam; we should consider if we can
// refactor into a single function that can produce both.
// In particular we note that bcf_open and bam_open are both
// just typedefs of hts_open. However, that could change.
SEXP load_bcf(SEXP bcf_file_r, SEXP index_file_r){
  if(TYPEOF(bcf_file_r) != STRSXP || length(bcf_file_r) != 1)
    error("bcf_file_r should contain a single string giving the bam file name");
  if(TYPEOF(index_file_r) != STRSXP || length(index_file_r) != 1)
    error("index_file_r should contain a single string giving the name of the bam index file");
  const char *bcf_file = CHAR(STRING_ELT( bcf_file_r, 0 ));
  const char *index_file = CHAR(STRING_ELT( index_file_r, 0 ));
  // make a bam_ptrs and load files and indices before returning these..
  // note that a vcfFile is a htsFile, but there bcfFile is not defined.
  vcfFile *vcf = bcf_open( bcf_file, "r" );
  if(!vcf)
    error("Unable to open specified bcf / vcf file");
  hts_idx_t *index = bcf_index_load2(bcf_file, index_file);
  if(!index){
    warning("No index file set; operations requiring index will not be available");
  }
  bcf_hdr_t *header = bcf_hdr_read(vcf);
  ////////// temporary block ////////////////////
  // look at the header structure a bit .. //////
  // This is just 0, 1, and 2.
  // We don't need to do this; we can use, "bcf_get_fmt( ... )
  // for each line. That may be a bit slower than pre-defining an array.
  /* int data_types[] = {BCF_DT_ID, BCF_DT_CTG, BCF_DT_SAMPLE}; */
  /* for(int i=0; i < 3; ++i){ */
  /*   for(int j=0; j < header->n[ data_types[i] ]; ++j){ */
  /*     Rprintf("\t%d : %s", j, header->id[data_types[i]][j].key); */
  /*   } */
  /*   Rprintf("\n"); */
  /* } */
  ///////////////////////////////////////////
  /// end of temp block.... ////////////////
  struct bcf_ptrs *ptr = malloc( sizeof(struct bcf_ptrs) );
  ptr->vcf = vcf;
  ptr->header = header;
  ptr->index = index;
  ptr->b_itr = 0;  // to set an iterator we have to make a range query.. 
  // create an external pointer with suitable tag and protect?
  SEXP tag = PROTECT( allocVector(STRSXP, 1) );
  SET_STRING_ELT( tag, 0, mkChar("bcf_ptrs") );
  // I'm not quite suire of the prot field, but might be useful
  SEXP prot = PROTECT( allocVector(STRSXP, 2) );
  SET_STRING_ELT( prot, 0, STRING_ELT(bcf_file_r, 0));
  SET_STRING_ELT( prot, 1, STRING_ELT(index_file_r, 0));
  SEXP ptr_r = PROTECT( R_MakeExternalPtr(ptr, tag, prot) );
  R_RegisterCFinalizerEx(ptr_r, finalise_bcf_ptrs, TRUE);
  UNPROTECT(3);
  return(ptr_r);
}

struct bam_ptrs* extract_bam_ptr(SEXP bam_ptr_r){
  if(TYPEOF(bam_ptr_r) != EXTPTRSXP)
    error("bam_ptr_r should be an external pointer");
  SEXP bam_tag = PROTECT(R_ExternalPtrTag(bam_ptr_r));
  // check the tag;
  if(TYPEOF(bam_tag) != STRSXP || length(bam_tag) != 1 || strcmp( CHAR(STRING_ELT(bam_tag, 0)), "bam_ptrs")){
    UNPROTECT(1);
    error("External pointer has incorrect tag");
    //    return(0);
  }
  struct bam_ptrs *bam = (struct bam_ptrs*)R_ExternalPtrAddr(bam_ptr_r);
  UNPROTECT(1);
  return(bam); // which should be checked by the caller
}

struct bcf_ptrs* extract_bcf_ptr(SEXP bcf_ptr_r){
  if(TYPEOF(bcf_ptr_r) != EXTPTRSXP)
    error("bcf_ptr_r should be an external pointer");
  SEXP bcf_tag = PROTECT(R_ExternalPtrTag(bcf_ptr_r));
  // check the tag;
  if(TYPEOF(bcf_tag) != STRSXP || length(bcf_tag) != 1 || strcmp( CHAR(STRING_ELT(bcf_tag, 0)), "bcf_ptrs")){
    UNPROTECT(1);
    error("External pointer has incorrect tag");
  }
  struct bcf_ptrs *bcf = (struct bcf_ptrs*)R_ExternalPtrAddr(bcf_ptr_r);
  UNPROTECT(1);
  return(bcf); // which should be checked by the caller
}

// returns 0 on error; Does not call error directly.
void *extract_ptr(SEXP ptr_r, const char* tag){
  if(TYPEOF(ptr_r) != EXTPTRSXP){
    warning("attempt to extract pointer from non pointer type");
    return(0);
  }
  SEXP tag_r = PROTECT(R_ExternalPtrTag(ptr_r));
  // check the tag;
  if(TYPEOF(tag_r) != STRSXP || length(tag_r) != 1 || strcmp( CHAR(STRING_ELT(tag_r, 0)), tag)){
    UNPROTECT(1);
    warning("External pointer has incorrect tag");
    return(0);
  }
  void *ptr = R_ExternalPtrAddr(ptr_r);
  UNPROTECT(1);
  return(ptr); // which should be checked by the caller
}


// Set an iterator for a specific region of a reference genome
// bam_ptr_r is an external pointer to a bam_ptrs pointer
// region_r is the name of a chromosome
// region_range_r is two integers
SEXP set_iterator(SEXP bam_ptr_r, SEXP region_r, SEXP region_range_r){
    struct bam_ptrs *bam = extract_bam_ptr(bam_ptr_r);
  if(!bam){
    // we could be clever here and try to recreate from the protect information
    // but we should think carefully about this
    error("External pointer is NULL");
  }
  if(!bam->sam || !bam->header || !bam->index)
    error("External pointer must contain header, sam and an index");
  if(TYPEOF(region_r) != STRSXP || length(region_r) < 1)
    error("region_r should be a character vector of positive length");
  if(TYPEOF(region_range_r) != INTSXP || length(region_range_r) != 2)
    error("region_range should be an integer vector of length 2");
  int *region_range = INTEGER(region_range_r);
  if(region_range[1] < region_range[0] || region_range[0] < 0)
    error("Unsorted region_range");
  const char *region = CHAR(STRING_ELT(region_r, 0));
  int target_id = sam_hdr_name2tid(bam->header, region);
  if(target_id < 0)
    error("Unable to find target id for specified region: %s", region);
  // if b_itr is not 0 then we should destroy it and set it to 0 before
  // we assign a new region.. (I don't actually know that this is needed,
  // but it seems a reasonable bet).
  if(bam->b_itr){
    hts_itr_destroy(bam->b_itr);
    bam->b_itr = 0;
  }
  bam->b_itr = sam_itr_queryi(bam->index, target_id, region_range[0], region_range[1]);
  return(bam_ptr_r);
}

SEXP clear_iterator(SEXP bam_ptr_r){
  struct bam_ptrs *bam = extract_bam_ptr(bam_ptr_r);
  if(!bam){
    // we could be clever here and try to recreate from the protect information
    // but we should think carefully about this
    error("External pointer is NULL");
  }
  if(bam->b_itr){
    hts_itr_destroy(bam->b_itr);
    bam->b_itr = 0;
  }
  return(bam_ptr_r);
}

// use iterator if not null
int read_next_entry( struct bam_ptrs *bam, bam1_t *b){
  if(bam->b_itr)
    return( sam_itr_next(bam->sam, bam->b_itr, b) );
  return( sam_read1(bam->sam, bam->header, b) );
}

// if( opt_flag & 1 ) ==> parse cigar data to give the ends of alignments
// min_mq_r: minimum mapping quality; only unique maximal mapping alignments
//           get non-0 scores. Setting this value to > 0 will remove such alignments.
//           This will make downstream analyses quicker, but removes some information.
SEXP build_query_index(SEXP bam_ptr_r, SEXP region_r, SEXP opt_flag_r, SEXP min_mq_r){
  struct bam_ptrs *bam = extract_bam_ptr(bam_ptr_r);
  if(!bam){
    // we could be clever here and try to recreate from the protect information
    // but we should think carefully about this
    error("External pointer is NULL");
  }
  if(TYPEOF(region_r) != STRSXP || length(region_r) < 1)
    error("region_r should be a character vector of positive length");
  if(TYPEOF(opt_flag_r) != INTSXP || length(opt_flag_r) != 1)
    error("opt_flag_r should be an integer vector containing one integer");
  int opt_flag = INTEGER(opt_flag_r)[0];

  if(TYPEOF(min_mq_r) != INTSXP || length(min_mq_r) != 1)
    error("min_mq_r should be an integer vector containing one integer");
  int min_mapq = INTEGER(min_mq_r)[0];
  
  if(!bam->index){
    error("External pointer does not contain an index structure");
  }
  const char *region = CHAR(STRING_ELT(region_r, 0));
  int target_id = sam_hdr_name2tid(bam->header, region);
  if(target_id < 0){
    warning("Unable to find target id for specified region: %s", region);
    return(R_NilValue);
  }
  int target_length = target_id < bam->header->n_targets ? bam->header->target_len[target_id] : 0;
  if(!target_length){
    warning("0 target length for specified region: %s", region);
    return(R_NilValue);
  }
  if(bam->qname_hash == 0){
    bam->qname_hash = init_sam_hash();
  }
  if(bam->qname_hash == 0)
    error("Unable to initialise a query name hash");

  // also check if we have the indexed_targets hash:
  if(bam->indexed_targets == 0)
    bam->indexed_targets = kh_init(int);

  if(bam->indexed_targets == 0)
    error("Unable to initialse a hash of indexed targets");

  int khash_ret;
  kh_put(int, bam->indexed_targets, target_id, &khash_ret);
  if(khash_ret == -1){
    warning("Failed to insert target id into bam->indexed_targets for: %d", target_id);
    return(R_NilValue);
  }
  if(khash_ret == 0){
    warning("Target %d already indexed. returning", target_id);
    return(R_NilValue);
  }
  
  hts_itr_t *b_itr = sam_itr_queryi(bam->index, target_id, 0, target_length);
  bam1_t *al = bam_init1();
  int ret;
  int count = 0;
  while((ret=sam_itr_next(bam->sam, b_itr, al)) >= 0){
    int flag = al->core.flag;
    // check that the segment is not unmapped; somehow it seems that unmapped segments
    // can be assigned locations;
    if(flag & 0x4)
      continue;
    int mapq = (int)al->core.qual;
    if(mapq < min_mapq)
      continue;
    int rpos_0 = 1 + (int)al->core.pos;
    int rpos_1 = rpos_0;
    int qpos_0 = 1;
    int qpos_1 = 1;
    int q_length = al->core.l_qseq;
    int qc_length = (opt_flag & 1) == 1 ? q_length : 0; 
    if((opt_flag & 1) && (al->core.n_cigar > 0)){
      uint32_t *cigar = bam_get_cigar(al);
      // 5 XOR [4, 5] --> 0, 1; this tests for BAM_CSOFT_CLIP and BAM_CHARD_CLIP
      qpos_0 += (bam_cigar_op(cigar[0]) ^ 5) < 2 ? bam_cigar_oplen( cigar[0] ) : 0;
      // If the first cigar op is hard clip, then it _should_ guarantee the presence
      // of a second cigar op; hence we should not need to check the cigar length here.
      qpos_0 += (bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP && bam_cigar_op(cigar[1]) == BAM_CSOFT_CLIP) ?
	bam_cigar_oplen(cigar[1]) : 0;
      qpos_1 = qpos_0;
      int op, op_length, op_type;
      // for this purpose the query should be incremented if the op is
      // M (0), I (1), =(7), X(8)
      // This can be tested by right shifting an integer that has those bits set
      // by the op value: i.e.
      // (0x0183 >> op) == 1 should check for query consumption in all cases. Thats a bit
      // more elegant than checking op_type and op for BAM_CHARD_CLIP.
      unsigned int query_bits = 0x0183;
      for(int i=0; i < al->core.n_cigar; ++i){
	op = bam_cigar_op(cigar[i]);
	op_length = bam_cigar_oplen(cigar[i]);
	op_type = bam_cigar_type(cigar[i]);
	rpos_1 += (op_type & CIG_OP_TP_RC) ? op_length : 0;
	qpos_1 += ((query_bits >> op) & 1) == 1 ? op_length : 0;
	//qpos_1 += (op_type == CIG_OP_TP_QC || op == BAM_CHARD_CLIP) ? op_length : 0;
	qc_length += (op == BAM_CHARD_CLIP) ? op_length : 0;
      }
    }
    srh_add_element(bam->qname_hash, bam_get_qname(al),
		    init_sam_record(target_id, rpos_0, rpos_1, flag,
				    qpos_0, qpos_1, q_length, mapq, qc_length));
    ++count;
  }
  SEXP ret_value = PROTECT(allocVector(INTSXP, 1));
  INTEGER(ret_value)[0] = count;
  UNPROTECT(1);
  return(ret_value);
}

SEXP query_positions(SEXP bam_ptr_r, SEXP query_ids_r){
  struct bam_ptrs *bam = extract_bam_ptr(bam_ptr_r);
  if(!bam){
    // we could be clever here and try to recreate from the protect information
    // but we should think carefully about this
    error("External pointer is NULL");
  }
  if(!bam->qname_hash)
    error("Query positions have not been indexed");

  if(TYPEOF(query_ids_r) != STRSXP || length(query_ids_r) < 1)
    error("query_ids_r should be a character vector of positive length");

  // lets try to make use of a kvec_t array here. It's already defined
  // in qname_hash.h
  sam_record_list query_pos;
  kv_init(query_pos);
  kvec_t(int) query_ids;
  kv_init(query_ids);
  for(int i=0; i < length(query_ids_r); ++i){
    sam_record_list* als = srh_get( bam->qname_hash, CHAR(STRING_ELT(query_ids_r, i)) );
    if(als){
      // kv_copy(sam_record, query_pos, *als);
      for(int j=0; j < als->n; ++j){
     	kv_push(sam_record, query_pos, kv_a(sam_record, *als, j));
	kv_push(int, query_ids, i+1);
      }
    } 
  }
  
  SEXP ret = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(ret, 0, allocVector(INTSXP, query_ids.n));
  SET_VECTOR_ELT(ret, 1, allocMatrix(INTSXP, SAM_RECORD_FIELDS_N, query_pos.n));
  SEXP row_names = mk_rownames( sam_record_field_names, SAM_RECORD_FIELDS_N );
  setAttrib(VECTOR_ELT(ret, 1), R_DimNamesSymbol, row_names);
  memcpy( INTEGER(VECTOR_ELT(ret, 0)), query_ids.a, sizeof(int) * query_ids.n );
  memcpy( INTEGER(VECTOR_ELT(ret, 1)), query_pos.a, sizeof(sam_record) * query_pos.n );
  kv_destroy(query_pos);
  kv_destroy(query_ids);
  UNPROTECT(1);
  return(ret);
}


// Parse a region in parallel using pthreads.
// Use 0-based counting for positions.
// The pointers to all data structures 
// initialised here. They should be 0 at this point.
void *alignments_region_thread(void *v_args){
  alignments_region_mt_args *args = (alignments_region_mt_args*)v_args;
  // First load the bam file and create the iterator.
  args->sam = sam_open(args->bam_file, "r");
  args->index = sam_index_load(args->sam, args->bam_index_file);
  if(!args->sam || !args->index)
    return(0);
  args->bam_it = sam_itr_queryi( args->index, args->target_id, args->start, args->end);
  
  uint32_t opt_flag = args->opts;
  size_t init_size = args->imatrix_initial_size;
  struct cigar_parse_options *cig_opts = &args->cig_opt;
  //  size_t region_length = (size_t)(args->end - args->start);
  // mandatory fields 
  args->al_core = init_i_matrix(AR_AL_RN, init_size);
  args->query_ids = init_str_array(init_size);
  args->al_ops = init_i_matrix(CIG_OPS_RN, args->imatrix_initial_size);
  // optional fields
  args->query_seq = init_str_array(init_size);
  args->query_qual = init_str_array(init_size);
  args->cigars = init_str_array(init_size);

  // optional fields held by args->cig_options: checking BOTH is
  // not a good idea; the information should be stored in one place only.
  if(bit_set(opt_flag, AR_Q_DEPTH) && cig_opts->depth)
    cig_opts->depth = cig_opts->depth + (args->start - cig_opts->depth_begin);

  if(bit_set(opt_flag, AR_Q_INTRON_DEPTH) && cig_opts->i_depth)
    cig_opts->i_depth = cig_opts->i_depth + (args->start - cig_opts->depth_begin);

  cig_opts->depth_l = args->end - args->start;
  cig_opts->depth_begin = args->start;

  // This is a bit wasteful because it creates the data structures even
  // if they are not needed. This is done because it simplifies what
  // needs to be done for check.
  args->mm_info = init_i_matrix(MM_INFO_RN, args->imatrix_initial_size);
  args->mate_core = init_i_matrix(AR_AL_MATE_RN, args->imatrix_initial_size);
  args->diff = init_i_matrix(Q_DIFF_RN, args->imatrix_initial_size);
  cig_opts->mm_info = bit_set(opt_flag, AR_AUX_MM) ? &args->mm_info : 0;
  cig_opts->diff = bit_set(opt_flag, AR_Q_DIFF) ? &args->diff : 0;
  
  // the iterator should be initialised before calling pthread_create
  int ret;
  bam1_t *al = bam_init1(); // needs to be destroyed at the end of parsing;
  int start = args->start;
  uint32_t *flag_filter = (uint32_t*)args->flag_filter;
  size_t seq_buffer_size=1024;
  // If we need to pass the query sequence to cigar_to_table; set up a buffer:
  int need_qseq = bit_set(opt_flag, AR_Q_SEQ) || bit_set(opt_flag, AR_Q_DIFF) || bit_set(opt_flag, AR_AUX_MM );
  cig_opts->qseq = need_qseq ? malloc(seq_buffer_size) : 0;
  struct cigar_string cigars_string = cigar_string_init(seq_buffer_size);
  int q_beg, q_end;
  size_t ops_begin_col, ops_end_col;
  while((ret=sam_itr_next(args->sam, args->bam_it, al)) >= 0){
    if(al->core.l_qseq < args->min_ql || al->core.qual < args->min_mq)
      continue;
    int rpos_0 = 1 + (int)al->core.pos;
    int do_push = rpos_0 >= start || cig_opts->include_left_als;
    int rpos_1 = rpos_0;
    uint32_t flag = al->core.flag;
    if((flag_filter[0] & flag) != flag_filter[0] || (flag_filter[1] & flag) > 0)
      continue;
    if(do_push){
      cig_opts->qseq_l = al->core.l_qseq;
      // bam_get_qual only returns a pointer to existing memory; hence it is cheap
      cig_opts->qqual = (char*)bam_get_qual(al); // This does not allocate memory

      // bam_seq_p, reads to seq_buffer (may need to change from const char* to char*)
      if(need_qseq)
	seq_buffer_size = bam_seq_p(al, &cig_opts->qseq, seq_buffer_size); 

      if(bit_set(opt_flag, AR_CIG))
	cigar_string_read( &cigars_string, al );
    }
    size_t al_i = 1 + args->al_core.col;
    // Parse the cigar data:
    int qcig_length = 0;
    cigar_to_table(al, al_i, &rpos_1, &args->al_ops,
		   &qcig_length, &q_beg, &q_end, &ops_begin_col, &ops_end_col, cig_opts, do_push);

    if(do_push){
      // Extract standard auxiliary fields;
      // it may be better to have a set of fields set by the user
      // that hasn't been done yet because it would require setting data structures
      // dynamically. (Not terribly difficult but).
      // Adding an option here to extract mq values if mate information has been requested
      // would be easy.
      int AS, NM, NH, IH;
      extract_int_aux_values(al, &AS, &NM, &NH, &IH);
      
      int exp_err = bit_set(opt_flag, AR_Q_ERR_EXP) && al->core.l_qseq && cig_opts->qqual ?
	estimate_read_errors((unsigned char*)cig_opts->qqual, al->core.l_qseq) : R_NaInt;
      // {0:flag, 1:r.beg, 2:r.end, 3:q.beg, 4:q.end, 5:mqual, 6:qlen, 7:qclen, 8:ops.0, 9:ops.1}
      push_column(&args->al_core, (int[AR_AL_RN]){flag, rpos_0, rpos_1, q_beg, q_end, (int)al->core.qual,
						    (int)al->core.l_qseq, qcig_length, 1+ops_begin_col,
						    ops_end_col, AS, NM, NH, IH, exp_err});
      str_array_push_cp(&args->query_ids, bam_get_qname(al));
      // optional fields:
      if(bit_set(opt_flag, AR_MT_INFO)) // add one to make R happy.
	push_column(&args->mate_core, (int[AR_AL_MATE_RN]){(int)al->core.mtid + 1, (int)al->core.mpos + 1, (int)al->core.isize} );
      if(bit_set(opt_flag, AR_Q_SEQ))
	str_array_push_cp(&args->query_seq, cig_opts->qseq);

      if(bit_set(opt_flag, AR_Q_QUAL))
	str_array_push_cp_n(&args->query_qual, cig_opts->qqual, al->core.l_qseq);

      if(bit_set(opt_flag, AR_CIG))
	str_array_push_cp(&args->cigars, cigars_string.buffer);
    }
  }
  // free up local resources
  free(cig_opts->qseq);
  cig_opts->qseq = 0;
  free(cigars_string.buffer);
  bam_destroy1(al);
  hts_idx_destroy(args->index);
  hts_itr_destroy(args->bam_it);
  sam_close(args->sam);
  pthread_exit((void*)args);
}

void *alignments_merge_thread(void *args_v){
  alignments_merge_args* args = (alignments_merge_args*)args_v;
  size_t int_size = sizeof(int);
  int *al_core_begin = args->al_core + (AR_AL_RN * args->core_off);
  int *al_ops_begin = args->al_ops + (CIG_OPS_RN * args->ops_off);
  int *diff_begin = args->diff ? args->diff + (Q_DIFF_RN * args->diff_off) : 0;
  int *mm_begin = args->mm ? args->mm + (MM_INFO_RN * args->mm_off) : 0;
  int *mate_begin = args->mate_core ? args->mate_core + (AR_AL_MATE_RN * args->mate_off) : 0;
  memcpy(al_core_begin,
	 args->t_args.al_core.data,
	 (AR_AL_RN * args->t_args.al_core.col * int_size));
  memcpy(al_ops_begin,
	 args->t_args.al_ops.data,
	 (CIG_OPS_RN * args->t_args.al_ops.col * int_size));
  if(args->diff)
    memcpy(diff_begin,
	   args->t_args.diff.data,
	   (Q_DIFF_RN * args->t_args.diff.col * int_size));
  if(args->mm)
    memcpy(mm_begin,
	   args->t_args.mm_info.data,
	   (MM_INFO_RN * args->t_args.mm_info.col * int_size));
  if(args->mate_core)
    memcpy(mate_begin, args->t_args.mate_core.data,
	   (AR_AL_MATE_RN * args->t_args.mate_core.col * int_size));
  
  /* // we then need to adjust the values for al_i in al_ops, diff and mm */
  for(size_t i=0; i < args->t_args.al_core.col; ++i){
    *(al_core_begin + AR_AL_RN * i + AR_AL_OPS_ROW) += args->ops_off;
    *(al_core_begin + AR_AL_RN * i + 1 + AR_AL_OPS_ROW) += args->ops_off;
  }
  for(size_t i=0; i < args->t_args.al_ops.col; ++i)
    *(al_ops_begin + CIG_OPS_RN  * i) += args->core_off;
  if(args->diff){
    for(size_t i=0; i < args->t_args.diff.col; ++i)
      *(diff_begin + Q_DIFF_RN * i) += args->core_off;
  }
  if(args->mm){
    for(size_t i=0; i < args->t_args.mm_info.col; ++i)
      *(mm_begin + MM_INFO_RN * i) += args->core_off;
  }
  pthread_exit((void*)args);
}

SEXP alignments_region_mt(SEXP bam_file_r, SEXP index_file_r,
			  SEXP region_r, SEXP region_range_r,
			  SEXP flag_filter_r, 
			  SEXP opt_flag_r, SEXP ref_seq_r,
			  SEXP min_mq_r, SEXP min_ql_r, SEXP max_intron_l_r,
			  SEXP n_threads_r){
  if(TYPEOF(bam_file_r) != STRSXP || length(bam_file_r) != 1)
    error("bam_file_r should contain a single string giving the bam file name");
  if(TYPEOF(index_file_r) != STRSXP || length(index_file_r) != 1)
    error("index_file_r should contain a single string giving the name of the bam index file");
  const char *bam_file = CHAR(STRING_ELT( bam_file_r, 0 ));
  const char *index_file = CHAR(STRING_ELT( index_file_r, 0 ));

  if(TYPEOF(n_threads_r) != INTSXP || length(n_threads_r) != 1)
    error("n_threads_r should be an integer vector of length 1");
  int n_threads = INTEGER(n_threads_r)[0];
  if(n_threads < 0){
    warning("negative number of threads (%d) specified. nthreads set to 1", n_threads);
    n_threads = 1;
  }
  if(n_threads > MAX_THREADS){
    warning("requested number (%d) of threads exceeds max allowed: using %d", n_threads, MAX_THREADS);
    n_threads = MAX_THREADS;
  }
  //// From HERE to PARSE_ARGS_END the code is identical to alignments_region
  //// I should refactor this part into smaller functions;
  //// could have functions like, extract_matrix_i(SEXP, *nrow, *ncol, *data, min_nrow, min_ncol)
  if(TYPEOF(region_r) != STRSXP || length(region_r) < 1)
    error("region_r should be a character vector of positive length");

  if(TYPEOF(region_range_r) != INTSXP || length(region_range_r) != 2)
    error("region_range should be an integer vector of length 2");
  int *region_range = INTEGER(region_range_r);

  if(region_range[1] <= region_range[0] || region_range[0] < 0)
    error("Unsorted region_range");
  const char *region = CHAR(STRING_ELT(region_r, 0));

  // open files temporarily to get a target id:
  samFile *tmp_sam = sam_open( bam_file, "r" );
  if(!tmp_sam)
    error("alignments_region_mt: unable to open bam file: %s", bam_file);
  hts_idx_t *tmp_index = sam_index_load(tmp_sam, index_file);
  if(!tmp_index){
    sam_close(tmp_sam);
    hts_idx_destroy(tmp_index);
    error("alignments_region_mt: Unable to load index from %s", index_file);
  }
  // we will use the header further below; destroy index and other later
  // at the same point. Note that any errors here have the potential to cause problems.
  sam_hdr_t *sam_header = sam_hdr_read(tmp_sam);
  int target_id = sam_hdr_name2tid(sam_header, region);
  // then destroy all the temporary values.
  sam_hdr_destroy(sam_header);
  sam_header = 0;
  hts_idx_destroy(tmp_index);
  sam_close(tmp_sam);
  
  if(target_id < 0){
    warning("Unable to find target id for specified region: %s", region);
    return(R_NilValue);
  }
  if(TYPEOF(flag_filter_r) != INTSXP || length(flag_filter_r) != 2)
    error("flag_filter should be an integer vector of length 2");
  if(TYPEOF(opt_flag_r) != INTSXP || length(opt_flag_r) != 1)
    error("opt_flag should be a single integer value");
  int opt_flag = asInteger(opt_flag_r);
  if(opt_flag < 0)
    opt_flag = 0;
  // cig_opt will be passed to cigar_to_table()
  // the default values are all 0s
  struct cigar_parse_options cig_opt = init_cigar_parse_options();
  cig_opt.depth_begin = region_range[0];
  
  if(bit_set(opt_flag, AR_Q_DIFF) && ((TYPEOF(ref_seq_r) != STRSXP) || length(ref_seq_r) != 1))
    error("Sequence divergence requested but ref_seq is not a character vector of length 1 ");
  const char *ref_seq = (TYPEOF(ref_seq_r) == STRSXP) ? CHAR(STRING_ELT(ref_seq_r, 0)) : 0;
  // which returns an int (which we can change to a size_t)
  size_t ref_seq_l = (size_t)(ref_seq != 0) ? LENGTH(STRING_ELT(ref_seq_r, 0)) : 0;
  cig_opt.rseq = ref_seq;
  cig_opt.rseq_l = ref_seq_l;
  int *flag_filter_rp = INTEGER(flag_filter_r);
  // default values; no filtering
  uint32_t flag_filter[2] = {0, 0};
  // override if set by user.
  if(flag_filter_rp[0] >= 0)
    flag_filter[0] = (uint32_t)flag_filter_rp[0];
  if(flag_filter_rp[1] >= 0)
    flag_filter[1] = (uint32_t)flag_filter_rp[1];

  if(TYPEOF(min_mq_r) != INTSXP || length(min_mq_r) != 1)
    error("min_mq_r should be an integer vector of length 1");
  if(TYPEOF(min_ql_r) != INTSXP || length(min_ql_r) != 1)
    error("min_ql_r should be an integer vector of length 1");
  int min_mq = asInteger(min_mq_r);
  int min_ql = asInteger(min_ql_r);

  // Set up the return data structure as a named list
  SEXP ret_names_r = PROTECT( mk_strsxp(ar_return_fields, AR_R_FIELDS_N) );
  SEXP ret_data = PROTECT(allocVector(VECSXP, length(ret_names_r)));
  setAttrib( ret_data, R_NamesSymbol, ret_names_r );

  // If depth requested we can add this to the data structure immediately..
  int *seq_depth = 0;
  /// I think I should get rid of the 1 + here.
  int region_length = 1 + region_range[1] - region_range[0];
  if(bit_set(opt_flag, AR_Q_DEPTH) | bit_set(opt_flag, AR_Q_INTRON_DEPTH))
    cig_opt.depth_l = region_length;
  if(bit_set(opt_flag, AR_Q_DEPTH)){
    SET_VECTOR_ELT( ret_data, 6, allocVector(INTSXP, region_length) );
    seq_depth = INTEGER(VECTOR_ELT(ret_data, 6));
    memset( seq_depth, 0, sizeof(int) * region_length );
    cig_opt.depth = seq_depth;
  }
  // And similarly for the intron depth:
  cig_opt.max_intron_length = MAX_INTRON_L; // unless set by user
  if(TYPEOF(max_intron_l_r) == INTSXP && length(max_intron_l_r) == 1){
    cig_opt.max_intron_length = INTEGER(max_intron_l_r)[0];
  }else{
    warning("max_intron_length set to default value (%d)", cig_opt.max_intron_length);
  }
  int *intron_depth = 0;
  if(bit_set(opt_flag, AR_Q_INTRON_DEPTH)){
    SET_VECTOR_ELT( ret_data, 10, allocVector(INTSXP, region_length));
    intron_depth = INTEGER(VECTOR_ELT(ret_data, 10));
    memset( intron_depth, 0, sizeof(int) * region_length );
    cig_opt.i_depth = intron_depth;
  }
  // PARSE_ARGS_END: end of (almost) identical code (copy pasted evil)

  pthread_t *threads = malloc(sizeof(pthread_t) * n_threads);
  alignments_region_mt_args *t_args = malloc(sizeof(alignments_region_mt_args) * n_threads);
  alignments_region_mt_args args_p = init_ar_args();
  args_p.bam_file = bam_file;
  args_p.bam_index_file = index_file;
  args_p.target_id = target_id;
  args_p.flag_filter = flag_filter;
  args_p.cig_opt = cig_opt;
  args_p.opts = opt_flag;
  args_p.min_mq = min_mq;
  args_p.min_ql = min_ql;
  int sub_range = region_length / n_threads;
  int start = region_range[0];
  int rem = region_length % n_threads;
  // Looks like I will need to create a new samFile for each of these
  for(int i=0; i < n_threads; ++i){
    t_args[i] = args_p; // the prototype
    t_args[i].start = start;
    t_args[i].end = start + sub_range + (i == 0 ? rem : 0);
    start = t_args[i].end;
    t_args[i].cig_opt.include_left_als = (i == 0) ? 1 : 0;
    pthread_create( &threads[i], NULL, alignments_region_thread, (void*)&t_args[i]);
  }
  void *status;
  int n_core=0, n_mate=0, n_ops=0, n_diff=0, n_mm=0;
  for(int i=0; i < n_threads; ++i){
    pthread_join(threads[i], &status);
  // Then merge the data and clean up stuff..
    n_core += t_args[i].al_core.col;
    n_mate += t_args[i].mate_core.col;
    n_ops += t_args[i].al_ops.col;
    n_diff += t_args[i].diff.col;
    n_mm += t_args[i].mm_info.col;
  }
  // depth (6) and intron depth have already been set (10)
  // There should be a better way of doing this, but..
  // The idea is to use the bit defined as the ret_data index
  // so instead of (ret_data, 0, region_r), we would use AR_<option name>
  // but this would require either reordering or other extensive changes
  // which have the potential to cause other problems. 
  SET_VECTOR_ELT(ret_data, 0, region_r);
  SET_VECTOR_ELT(ret_data, 1, allocVector(STRSXP, n_core));
  SET_VECTOR_ELT(ret_data, 2, allocMatrix(INTSXP, AR_AL_RN, n_core));
  setAttrib(VECTOR_ELT(ret_data, 2), R_DimNamesSymbol, mk_rownames(ar_al_rownames, AR_AL_RN));
  SET_VECTOR_ELT(ret_data, 3, allocMatrix(INTSXP, CIG_OPS_RN, n_ops));
  setAttrib(VECTOR_ELT(ret_data, 3), R_DimNamesSymbol, mk_rownames(cig_ops_rownames, CIG_OPS_RN));
  if(bit_set(opt_flag, AR_Q_SEQ))
    SET_VECTOR_ELT(ret_data, 4, allocVector(STRSXP, n_core));
  if(bit_set(opt_flag, AR_Q_DIFF)){
    SET_VECTOR_ELT(ret_data, 5, allocMatrix(INTSXP, Q_DIFF_RN, n_diff));
    setAttrib(VECTOR_ELT(ret_data, 5), R_DimNamesSymbol, mk_rownames(q_diff_rownames, Q_DIFF_RN));
  }
  if(bit_set(opt_flag, AR_CIG))
    SET_VECTOR_ELT(ret_data, 7, allocVector(STRSXP, n_core));
  if(bit_set(opt_flag, AR_Q_QUAL))
    SET_VECTOR_ELT(ret_data, 8, allocVector(STRSXP, n_core));
  if(bit_set(opt_flag, AR_AUX_MM)){
    SET_VECTOR_ELT(ret_data, 9, allocMatrix(INTSXP, MM_INFO_RN, n_mm));
    setAttrib(VECTOR_ELT(ret_data, 9), R_DimNamesSymbol, mk_rownames(mm_info_rownames, MM_INFO_RN));
  }
  if(bit_set(opt_flag, AR_MT_INFO)){
    SET_VECTOR_ELT(ret_data, 11, allocMatrix(INTSXP, AR_AL_MATE_RN, n_mate));
    setAttrib(VECTOR_ELT(ret_data, 11), R_DimNamesSymbol, mk_rownames(al_mate_rownames, AR_AL_MATE_RN));
  }
  
  // The following (with lots of copying of strings), may be wortwhile to
  // also thread. It should be fairly straightforward.
  int *core_al = INTEGER(VECTOR_ELT(ret_data, 2));
  int *mate_core = bit_set(opt_flag, AR_MT_INFO) ? INTEGER(VECTOR_ELT(ret_data, 11)) : 0;
  int *cig_ops = INTEGER(VECTOR_ELT(ret_data, 3));
  int *diff = bit_set(opt_flag, AR_Q_DIFF) ? INTEGER(VECTOR_ELT(ret_data, 5)) : 0;
  int *aux_mm = bit_set(opt_flag, AR_AUX_MM) ? INTEGER(VECTOR_ELT(ret_data, 9)) : 0;

  SEXP query_id = PROTECT(VECTOR_ELT(ret_data, 1));
  SEXP query_seq = PROTECT(VECTOR_ELT(ret_data, 4));
  SEXP cigars = PROTECT(VECTOR_ELT(ret_data, 7));
  SEXP qqual = PROTECT(VECTOR_ELT(ret_data, 8));
  
  // First check the single threaded way:
  //  int core_i=0, ops_i=0, diff_i=0, mm_i=0;
  size_t core_off = 0, mate_off=0, ops_off=0, diff_off=0, mm_off=0;
  alignments_merge_args *merge_args = malloc(sizeof(alignments_merge_args) * n_threads);
  alignments_merge_args merge_args_p = init_ar_merge_args();
  merge_args_p.al_core = core_al;
  merge_args_p.mate_core = mate_core;
  merge_args_p.al_ops = cig_ops;
  merge_args_p.diff = diff;
  merge_args_p.mm = aux_mm;
  // the checks in the two following statements are redundant as the
  // value of diff and aux_mm should be 0 if the bit is not set.
  /* merge_args_p.diff = bit_set(opt_flag, AR_Q_DIFF) ? diff : 0; */
  /* merge_args_p.mm = bit_set(opt_flag, AR_AUX_MM) ? aux_mm : 0; */
  pthread_t *m_threads = malloc(sizeof(pthread_t) * n_threads);
  // set up the imatrix merge threads:
  for(int i=0; i < n_threads; ++i){
    merge_args[i] = merge_args_p;
    merge_args[i].t_args = t_args[i];
    arm_set_offsets(&merge_args[i], core_off, ops_off, diff_off, mm_off, mate_off);
    pthread_create( &m_threads[i], NULL, alignments_merge_thread, (void*)&merge_args[i] );
    core_off += t_args[i].al_core.col;
    ops_off +=  t_args[i].al_ops.col;
    diff_off += t_args[i].diff.col;
    mm_off += t_args[i].mm_info.col;
    mate_off += t_args[i].mate_core.col;
  }
  
  core_off = ops_off = diff_off = mm_off = 0;
  for(int i=0; i < n_threads; ++i){

    for(int j=0; j < t_args[i].al_core.col; ++j){
      SET_STRING_ELT(query_id, core_off + j, mkChar(t_args[i].query_ids.strings[j]));
      if(bit_set(opt_flag, AR_Q_SEQ))
	SET_STRING_ELT(query_seq, core_off + j, mkChar(t_args[i].query_seq.strings[j]));
      if(bit_set(opt_flag, AR_CIG))
	SET_STRING_ELT(cigars, core_off + j, mkChar(t_args[i].cigars.strings[j]));
      if(bit_set(opt_flag, AR_Q_QUAL))
	SET_STRING_ELT(qqual, core_off + j, mkChar(t_args[i].query_qual.strings[j]));
    }
    core_off += t_args[i].al_core.col;
    ops_off +=  t_args[i].al_ops.col;
    diff_off += t_args[i].diff.col;
    mm_off += t_args[i].mm_info.col;

    // and free data: // this now done by free_ar_args further down
    /* str_array_free( &t_args[i].query_ids ); */
    /* if(bit_set(opt_flag, AR_Q_SEQ)) */
    /*   str_array_free( &t_args[i].query_seq ); */
    /* if(bit_set(opt_flag, AR_CIG)) */
    /*   str_array_free( &t_args[i].cigars ); */
    /* if(bit_set(opt_flag, AR_Q_QUAL)) */
    /*   str_array_free( &t_args[i].query_qual ); */
  }
  // Join the other threads
  for(int i=0; i < n_threads; ++i){
    pthread_join(m_threads[i], &status);
    free_ar_args( &t_args[i] );
  }
  /* for(int i=0; i < n_threads; ++i){ */
  /*   // free data from t_args[i] */
  /*   free( t_args[i].al_core.data ); */
  /*   free( t_args[i].al_ops.data ); */
  /*   if(bit_set(opt_flag, AR_AUX_MM)) */
  /*     free( t_args[i].mm_info.data ); */
  /*   if(bit_set(opt_flag, AR_Q_DIFF)) */
  /*     free( t_args[i].diff.data ); */
  /* } */
  free(t_args);
  free(merge_args);
  free(threads);
  free(m_threads);
  
  UNPROTECT(6);
  return(ret_data);
}

// region_r: the region as:  ref_name:beg-end
//                           or just ref
// region_range_r: integer vector with start and end coordinates
// range_r: begin and end locations on the region
// bam_ptr_r: an external pointer with tag bam_ptrs and a pointer
//            to a bam_ptrs struct
// flag_filter_r: Two integer values:
//            1. need: alignments must have all bits set in flag
//            2. exclude: alignments must have none of the bits set in flag
//  A negative value is treated as a default value:
//            1. 0
//            2. 0
//  opt_flag_r
//  Bitwise OR Flag giving options about what data to return.. 
//     bit 1 = seq_data as a character vector array.
//     bit 2 = return positions that differ from reference. Requires that reference is given
//     bit 3 = calculate the coverage depth
//     bit 4 = reconstruct a cigar string
//     bit 5 = return base qualities
//     bit 6 = return mate information
//     bit 7 = parse MM auxiliary string.
//     bit 8 = return intron depth
//  min_mq_r : A minimum mapping quality
//  min_ql_r : A minimum query length, taken from bam1_core_t.qlen
//             which does not necessarily equate to the read length
//  max_intron_l : max intron length. Only relevant if intron depth is requested
SEXP alignments_region(SEXP region_r, SEXP region_range_r,
		       SEXP bam_ptr_r, SEXP flag_filter_r, 
		       SEXP opt_flag_r, SEXP ref_seq_r,
		       SEXP min_mq_r, SEXP min_ql_r, SEXP max_intron_l_r){
  if(TYPEOF(region_r) != STRSXP || length(region_r) < 1)
    error("region_r should be a character vector of positive length");
  struct bam_ptrs *bam = extract_bam_ptr(bam_ptr_r);
  if(!bam){
    // we could be clever here and try to recreate from the protect information
    // but we should think carefully about this
    error("External pointer is NULL");
  }
  if(!bam->index){
    error("External pointer does not contain an index structure");
  }
  if(TYPEOF(region_range_r) != INTSXP || length(region_range_r) != 2)
    error("region_range should be an integer vector of length 2");
  int *region_range = INTEGER(region_range_r);
  if(region_range[1] <= region_range[0] || region_range[0] < 0)
    error("Unsorted region_range");
  const char *region = CHAR(STRING_ELT(region_r, 0));
  int target_id = sam_hdr_name2tid(bam->header, region);
  if(target_id < 0){
    warning("Unable to find target id for specified region: %s", region);
    return(R_NilValue);
  }
  if(TYPEOF(flag_filter_r) != INTSXP || length(flag_filter_r) != 2)
    error("flag_filter should be an integer vector of length 2");
  if(TYPEOF(opt_flag_r) != INTSXP || length(opt_flag_r) != 1)
    error("opt_flag should be a single integer value");
  int opt_flag = asInteger(opt_flag_r);
  if(opt_flag < 0)
    opt_flag = 0;
  // cig_opt will be passed to cigar_to_table()
  // the default values are all 0s
  struct cigar_parse_options cig_opt = init_cigar_parse_options();
  cig_opt.depth_begin = region_range[0];
  
  if(bit_set(opt_flag, AR_Q_DIFF) && ((TYPEOF(ref_seq_r) != STRSXP) || length(ref_seq_r) != 1))
    error("Sequence divergence requested but ref_seq is not a character vector of length 1 ");
  const char *ref_seq = (TYPEOF(ref_seq_r) == STRSXP) ? CHAR(STRING_ELT(ref_seq_r, 0)) : 0;
  // which returns an int (which we can change to a size_t)
  size_t ref_seq_l = (size_t)(ref_seq != 0) ? LENGTH(STRING_ELT(ref_seq_r, 0)) : 0;
  cig_opt.rseq = ref_seq;
  cig_opt.rseq_l = ref_seq_l;
  int *flag_filter_rp = INTEGER(flag_filter_r);
  // default values; no filtering
  uint32_t flag_filter[2] = {0, 0};
  // override if set by user.
  if(flag_filter_rp[0] >= 0)
    flag_filter[0] = (uint32_t)flag_filter_rp[0];
  if(flag_filter_rp[1] >= 0)
    flag_filter[1] = (uint32_t)flag_filter_rp[1];

  if(TYPEOF(min_mq_r) != INTSXP || length(min_mq_r) != 1)
    error("min_mq_r should be an integer vector of length 1");
  if(TYPEOF(min_ql_r) != INTSXP || length(min_ql_r) != 1)
    error("min_ql_r should be an integer vector of length 1");
  int min_mq = asInteger(min_mq_r);
  int min_ql = asInteger(min_ql_r);

  // Set up the return data structure as a named list
  SEXP ret_names_r = PROTECT( mk_strsxp(ar_return_fields, AR_R_FIELDS_N) );
  SEXP ret_data = PROTECT(allocVector(VECSXP, length(ret_names_r)));
  setAttrib( ret_data, R_NamesSymbol, ret_names_r );

  // If depth requested we can add this to the data structure immediately..
  int *seq_depth = 0;
  int region_length = 1 + region_range[1] - region_range[0];
  if(bit_set(opt_flag, AR_Q_DEPTH)){
    SET_VECTOR_ELT( ret_data, 6, allocVector(INTSXP, region_length) );
    seq_depth = INTEGER(VECTOR_ELT(ret_data, 6));
    memset( seq_depth, 0, sizeof(int) * region_length );
    cig_opt.depth = seq_depth;
    cig_opt.depth_l = region_length;
  }
  // And similarly for the intron depth:
  int max_intron_length = MAX_INTRON_L; // unless set by user
  if(TYPEOF(max_intron_l_r) == INTSXP && length(max_intron_l_r) == 1){
    max_intron_length = INTEGER(max_intron_l_r)[0];
  }else{
    warning("max_intron_length set to default value (%d)", max_intron_length);
  }
  int *intron_depth = 0;
  cig_opt.max_intron_length = max_intron_length;
  if(bit_set(opt_flag, AR_Q_INTRON_DEPTH)){
    SET_VECTOR_ELT( ret_data, 10, allocVector(INTSXP, region_length));
    intron_depth = INTEGER(VECTOR_ELT(ret_data, 10));
    memset( intron_depth, 0, sizeof(int) * region_length );
    cig_opt.i_depth = intron_depth;
    cig_opt.depth_l = region_length;
  }
  // We use sam_itr_queryi, so that we can specify the region
  // numerically thus also allowing the same function to return the sequence depth for the region.
  hts_itr_t *b_itr = sam_itr_queryi(bam->index, target_id, region_range[0], region_range[1]);
  // al should be destroyed with bam_destroy1
  bam1_t *al = bam_init1();
  int r=0;
  int al_count = 0;
  // Data structures that will hold the return data:
  size_t init_size = OPS_INIT_SIZE;
  struct str_array query_ids = init_str_array(init_size);
  // query_seq is filled ony if opt_flag & (1 << AR_Q_SEQ) == 1
  struct str_array query_seq = init_str_array(init_size);
  // query_qual is filled if AR_Q_QUAL is set
  struct str_array query_qual = init_str_array(init_size);
  // We unfortunaly have ended up with a struct cigar_string,
  // and cigar, cigars and cigars_string variables; That is BAD, and ought to be cleaned up
  struct str_array cigars = init_str_array(init_size);
  struct cigar_string cigars_string = cigar_string_init(init_size);
  
  struct i_matrix al_coord = init_i_matrix( CIG_OPS_RN, init_size );
  SEXP al_coord_rownames_r = PROTECT(mk_rownames( cig_ops_rownames, CIG_OPS_RN ));

  // For sequence variants in the reads:
  //  struct i_matrix var_coord = init_i_matrix( 4, init_size );
  //  SEXP var_coord_row_names_r = PROTECT(mk_rownames( (const char*[]){"al.i", "r.pos", "q.pos", "nuc"}, 4 ));
  struct i_matrix var_coord = init_i_matrix( Q_DIFF_RN, init_size );
  SEXP var_coord_row_names_r = PROTECT(mk_rownames( q_diff_rownames, Q_DIFF_RN ));
  if(bit_set(opt_flag, AR_Q_DIFF))
    cig_opt.diff = &var_coord;
  // For base modification data:
  struct i_matrix mm_info = init_i_matrix( MM_INFO_RN, init_size );
  SEXP mm_info_rownames_r = PROTECT(mk_rownames( mm_info_rownames, MM_INFO_RN ));
  if(bit_set(opt_flag, AR_AUX_MM))
    cig_opt.mm_info = &mm_info;
  
  // A structure holding: flag, r.beg, r.end, q.beg, q.end, mqual, qlen, qclen, ops.0, ops.1
  struct i_matrix al_vars = init_i_matrix( AR_AL_RN, init_size );
  SEXP al_var_dimnames = PROTECT( mk_rownames(ar_al_rownames, AR_AL_RN));

  // qqual and qseq are a bit weird; memory may be allocated, but if so, the pointers will be
  // held in a str_array object and will be deallocated after these have been copied to the
  // return data structure. 
  char *qseq = 0; 
  // Like qseq, but for the quality values. These can be obtained using
  // bam_get_qual() which returns a pointer to the qual data
  //                of length bam->core.l_qseq
  unsigned char *qqual=0; //
  
  while((r=sam_itr_next(bam->sam, b_itr, al)) >= 0){
    qseq = 0, qqual = 0;
    uint32_t flag = (uint32_t)al->core.flag;
    // any bits set in flag_filter[0] must also be set in flag
    if((flag & flag_filter[0]) != flag_filter[0])
      continue;
    // any bits set in flag_filter[1] must be 0 in flag
    if((flag & flag_filter[1]) > 0)
      continue;
    if(al->core.qual < min_mq)
      continue;
    if(al->core.l_qseq < min_ql)
      continue;
    str_array_push_cp(&query_ids, bam_get_qname(al));
    // read the cigar and push it to cigars
    if(bit_set(opt_flag, AR_CIG)){
      cigar_string_read( &cigars_string, al );
      str_array_push_cp(&cigars, cigars_string.buffer );
    }
    if(bit_set(opt_flag, AR_Q_SEQ) || bit_set(opt_flag, AR_Q_DIFF)){
      qseq = bam_seq(al);
      qqual = bam_get_qual(al);
      if(bit_set(opt_flag, AR_Q_SEQ)){
	if(qseq)
	  str_array_push(&query_seq, qseq);
	else
	  str_array_push_cp(&query_seq, "*");
      }
    }
    // POTENTIAL PROBLEM; the quality is stored as a 0 terminated sequence
    // of bytes. But in theory it can contain 0 values; To avoid this I would
    // have to increment every value by a set amount (eg. 33); But this is slow.
    // I will need some time to consider this.
    // In R I don't think I have a way to handle the data encoded in this way.
    if(bit_set(opt_flag, AR_Q_QUAL)){
      if(qqual == 0)
	qqual = bam_get_qual(al);
      if(qqual)
	str_array_push_cp_n(&query_qual, (const char*)qqual, al->core.l_qseq);
      else
	str_array_push_cp(&query_qual, "*");
    }
    al_count++;

    // then go through the cigar string and for each op create a new column
    // of offsets for the al_count
    //    uint32_t *cigar = bam_get_cigar(al);
    //    int32_t cigar_l = al->core.n_cigar;
    int r_pos = 1 + (int)al->core.pos;
    int r_pos_0 = r_pos;
    //    int q_pos = 1;
    int q_beg=-1;
    int q_end=0;
    int qcig_length = 0;
    size_t ops_begin_col;
    size_t ops_end_col;
    cig_opt.qseq = qseq;
    cig_opt.qqual = (char*)qqual;
    cig_opt.qseq_l = qseq ? al->core.l_qseq : 0;
    // Parse cigar data to R tabular data
    cigar_to_table(al, al_count, &r_pos, &al_coord,
		   &qcig_length, &q_beg, &q_end, &ops_begin_col, &ops_end_col, &cig_opt, 1);

    // Extract standard auxiliary fields; note this is still experimental;
    int AS, NM, NH, IH;
    extract_int_aux_values(al, &AS, &NM, &NH, &IH);

    int exp_err = bit_set(opt_flag, AR_Q_ERR_EXP) && al->core.l_qseq && qqual ?
      estimate_read_errors(qqual, al->core.l_qseq) : R_NaInt;
    
    push_column(&al_vars, (int[AR_AL_RN]){(int)al->core.flag, r_pos_0, r_pos, q_beg, q_end, (int)al->core.qual,
					    (int)al->core.l_qseq, qcig_length, 1+(int)ops_begin_col, (int)ops_end_col,
					    AS, NM, NH, IH, exp_err});
    // We should modify the function so that we do not malloc and free for every query
    // sequence.
    // qseq may be 0, but that should be safe to free.
    // We only need to free qseq if we did not store it in the str_array_object
    // In fact we should never need to do this.
    if(!bit_set(opt_flag, AR_Q_SEQ) && qseq )
      free(qseq);
    qseq = 0;
  }
  // Destroy the iterator as we will have no further use of it:
  hts_itr_destroy( b_itr );
  bam_destroy1( al );
  // query ids
  SET_VECTOR_ELT(ret_data, 0, region_r);
  SET_VECTOR_ELT(ret_data, 1, allocVector(STRSXP, query_ids.length));
  if(bit_set(opt_flag, AR_CIG))
    SET_VECTOR_ELT(ret_data, 7, allocVector(STRSXP, cigars.length));
  // Handle query ids and query sequences at the same time as they have
  // the same length.
  if(bit_set(opt_flag, AR_Q_SEQ))
    SET_VECTOR_ELT(ret_data, 4, allocVector(STRSXP, query_ids.length));
  if(bit_set(opt_flag, AR_Q_QUAL))
    SET_VECTOR_ELT(ret_data, 8, allocVector(STRSXP, query_ids.length));
  SEXP q_ids_r = VECTOR_ELT( ret_data, 1 );
  SEXP cigars_r = VECTOR_ELT( ret_data, 7 );
  for(size_t i=0; i < query_ids.length; ++i){
    SET_STRING_ELT( q_ids_r, i, mkChar( query_ids.strings[i] ));
    if(bit_set(opt_flag, AR_CIG))
      SET_STRING_ELT( cigars_r, i, mkChar( cigars.strings[i] ));
    if(bit_set(opt_flag, AR_Q_SEQ))
      SET_STRING_ELT( VECTOR_ELT(ret_data, 4), i, mkChar( query_seq.strings[i] ));
    if(bit_set(opt_flag, AR_Q_QUAL))
      SET_STRING_ELT( VECTOR_ELT(ret_data, 8), i, mkChar( query_qual.strings[i] ));
  }
  str_array_free(&query_ids);
  str_array_free(&query_seq);
  str_array_free(&query_qual);
  str_array_free(&cigars);
  cigar_string_free(&cigars_string);
  // alignment variables
  SET_VECTOR_ELT(ret_data, 2, allocMatrix(INTSXP, al_vars.nrow, al_vars.col));
  setAttrib( VECTOR_ELT(ret_data, 2), R_DimNamesSymbol, al_var_dimnames );
  memcpy( INTEGER(VECTOR_ELT(ret_data, 2)), al_vars.data, sizeof(int) * al_vars.nrow * al_vars.col );
  free(al_vars.data);
  // alignment coordinates
  SET_VECTOR_ELT(ret_data, 3, allocMatrix(INTSXP, al_coord.nrow, al_coord.col));
  //  setAttrib( VECTOR_ELT(ret_data, 1), R_RowNamesSymbol, rowNames_r);
  setAttrib( VECTOR_ELT(ret_data, 3), R_DimNamesSymbol, al_coord_rownames_r);
  memcpy( INTEGER(VECTOR_ELT(ret_data, 3)), al_coord.data, sizeof(int) * al_coord.nrow * al_coord.col);
  free(al_coord.data);
  // If we have variant data then we need to do something with it before freeing up the used memory.
  if(bit_set(opt_flag, AR_Q_DIFF)){
    SET_VECTOR_ELT(ret_data, 5, allocMatrix(INTSXP, var_coord.nrow, var_coord.col));
    setAttrib( VECTOR_ELT(ret_data, 5), R_DimNamesSymbol, var_coord_row_names_r );
    memcpy( INTEGER(VECTOR_ELT(ret_data, 5)), var_coord.data, sizeof(int) * var_coord.nrow * var_coord.col );
  }
  if(bit_set(opt_flag, AR_AUX_MM)){
    SET_VECTOR_ELT(ret_data, 9, allocMatrix(INTSXP, mm_info.nrow, mm_info.col));
    setAttrib( VECTOR_ELT(ret_data, 9), R_DimNamesSymbol, mm_info_rownames_r);
    memcpy( INTEGER(VECTOR_ELT(ret_data, 9)), mm_info.data, sizeof(int) * mm_info.nrow * mm_info.col);
  }
  free(var_coord.data);
  free(mm_info.data);
  UNPROTECT(6);
  return(ret_data);
}

SEXP count_region_alignments(SEXP region_r, SEXP region_range_r,
			     SEXP bam_ptr_r, SEXP flag_filter_r, 
			     SEXP min_mq_r){
  if(TYPEOF(region_r) != STRSXP || length(region_r) < 1)
    error("region_r should be a character vector of positive length");
  struct bam_ptrs *bam = extract_bam_ptr(bam_ptr_r);
  if(!bam){
    // we could be clever here and try to recreate from the protect information
    // but we should think carefully about this
    error("External pointer is NULL");
  }
  if(!bam->index){
    error("External pointer does not contain an index structure");
  }
  if(TYPEOF(region_range_r) != INTSXP || length(region_range_r) != 2)
    error("region_range should be an integer vector of length 2");
  int *region_range = INTEGER(region_range_r);
  if(region_range[1] < region_range[0] || region_range[0] < 0)
    error("Unsorted region_range");
  const char *region = CHAR(STRING_ELT(region_r, 0));
  int target_id = sam_hdr_name2tid(bam->header, region);
  if(target_id < 0){
    warning("Unable to find target id for specified region: %s", region);
    return(R_NilValue);
  }
  if(TYPEOF(min_mq_r) != INTSXP || length(min_mq_r) != 1)
    error("min_mq_r should be an integer vector of length 1");
  int min_mq = asInteger(min_mq_r);
  //  int min_ql = asInteger(min_ql_r);

  if(TYPEOF(flag_filter_r) != INTSXP || length(flag_filter_r) != 2)
    error("flag_filter should be an integer vector of length 2");
  int *flag_filter_rp = INTEGER(flag_filter_r);
  // default values; no filtering
  uint32_t flag_filter[2] = {0, 0};
  // override if set by user.
  if(flag_filter_rp[0] >= 0)
    flag_filter[0] = (uint32_t)flag_filter_rp[0];
  if(flag_filter_rp[1] >= 0)
    flag_filter[1] = (uint32_t)flag_filter_rp[1];

  // We use sam_itr_queryi, so that we can specify the region
  // numerically thus also allowing the same function to return the sequence depth for the region.
  hts_itr_t *b_itr = sam_itr_queryi(bam->index, target_id, region_range[0], region_range[1]);
  bam1_t *al = bam_init1();
  int r=0;
  int n=0; // the number of alignments
  while((r=sam_itr_next(bam->sam, b_itr, al)) >= 0){
    uint32_t flag = (uint32_t)al->core.flag;
    // any bits set in flag_filter[0] must also be set in flag
    if((flag & flag_filter[0]) != flag_filter[0])
      continue;
    // any bits set in flag_filter[1] must be 0 in flag
    if((flag & flag_filter[1]) > 0)
      continue;
    if(al->core.qual < min_mq)
      continue;
    ++n;
  }
  // It might be good to call hts_itr_destroy here to avoid a memory leak:
  // hts_itr_destroy( b_itr );
  SEXP ret_value = PROTECT(allocVector(INTSXP, 1));
  *(INTEGER(ret_value)) = n;
  UNPROTECT(1);
  return( ret_value );
}

// Reads from an unindexed sam / bam / cram file
// bam_ptr_r : an external pointer to a bam_ptr struct
// n_r       : the number of reads to read (a single integer)
// ret_f_r   : a bitwise flags giving the columns to return data from
//            (a single integer)
// sel_flags : flags which select which alignmens are returned
//            three integer values giving:
//            f : only include reads with all of the bits set in flag
//            F : exclude reads with any of the bits set
//            q : minimum mapping quality 
// This mimics some of the options in samtools view, but 
// omits, -r (by readgroup) and selection based on location in the file
SEXP sam_read_n(SEXP bam_ptr_r, SEXP n_r, SEXP ret_f_r,
		SEXP sel_flags_r){
  // extract_bam_ptr does the error checking and handling
  struct bam_ptrs *bam = extract_bam_ptr(bam_ptr_r);
  if(!bam || !bam->sam){
    // we could be clever here and try to recreate from the protect information
    // but we should think carefully about this
    error("External pointer or sam file is NULL");
  }
  if(TYPEOF(n_r) != INTSXP || length(n_r) != 1)
    error("n_r should be an integer vector of length 1");
  if(TYPEOF(ret_f_r) != INTSXP || length(ret_f_r) != 1)
    error("ret_f_r should be an integer vector of length 1");
  if(TYPEOF(sel_flags_r) != INTSXP || length(sel_flags_r) != 3)
    error("sel_flags_r should be an integer vector of length 3 (f, F, q)");
  int n=asInteger(n_r);
  unsigned int ret_flag = asInteger(ret_f_r);
  int *sel_flags = INTEGER(sel_flags_r);
  unsigned int f_flag = (unsigned int)sel_flags[0];
  unsigned int F_flag = (unsigned int)sel_flags[1];
  int min_q = sel_flags[2];
  if(n <= 0)
    error("you must request at least one entry");
  // The dimensions and rownames are defined in common.h
  SEXP ret_data_names = PROTECT(mk_strsxp( srn_return_fields, SRN_FIELDS_N ));
  SEXP ret_data = PROTECT( allocVector(VECSXP, length(ret_data_names)) );
  setAttrib( ret_data, R_NamesSymbol, ret_data_names );
  // The last element is a single integer counting the number of entries returned:
  SET_VECTOR_ELT(ret_data, length(ret_data_names)-1, allocVector(INTSXP, 1));
  unsigned int flag = 1;
  // elements up to S_AUX are simple vectors with type defined by R_ret_types
  for(int i=0; i <= S_AUX; ++i){
    if(ret_flag & flag)
      SET_VECTOR_ELT(ret_data, i, allocVector( R_ret_types[i], n ));
    flag = flag << 1;
  }
  int *count = INTEGER(VECTOR_ELT(ret_data, length(ret_data_names)-1));
  *count = 0;
  // Elements after S_AUX, may take other forms; in general they will be
  // matrices.
  // The last element itself is just a single integer giving the number of
  // alignments parsed.
  // These additional matrices have unknown sizes and hence have to be copied
  // into the appropriate R data structures. This is to avoid using the R SEXP
  // resize function (I forget the name here). That anyway doesn't do anything clever
  // to avoid unnecessary memory allocations.
  bam1_t *b = bam_init1();
  size_t seq_buffer_size = 500;
  char *seq_buffer = malloc(seq_buffer_size);
  kstring_t aux_str = KS_INITIALIZE;
  // If S_AUX_MM is set then we need to parse the cigar table. We don't have to return
  // it to the R session, but for now it is easier to simply force this behaviour.
  // Alternatively, if the cigar is not parsed don't set reference positions.
  // We can work out the best option later. But for now we will force the behaviour.
  if( ret_flag & (1 << S_AUX_MM) )
    ret_flag |= (1 << S_CIG_TABLE);
  // If S_CIG_TABLE is set, then initialise to i_matrix arrays;
  // one will hold q_info (the actual query lengths, beg and end positions before clipping)
  // and the other will hold the cigar operations in a table of unknown length;
  int ret_cig_ops = ret_flag & (1 << S_CIG_TABLE);
  struct i_matrix al_coord = init_i_matrix(CIG_OPS_RN, (ret_cig_ops ? OPS_INIT_SIZE : 0));
  // they query_info will hold the q_length, q_cig_length, q_beg, q_end. (possibly also
  struct i_matrix query_info = init_i_matrix(Q_INFO_RN, (ret_cig_ops ? n : 0));
  // The al_coord.data and query_info.data should be freed after copying to an R SEXP

  // If S_AUX_MM specified, then parse MM information into an i_matrix that will need to be copied.
  // the 6 rows of the i_matrix are: al_i, query_pos, mod_code, mod_n, mod_likelihood, ref_pos
  struct i_matrix mm_info = init_i_matrix(MM_INFO_RN, (ret_flag & (1 << S_AUX_MM)) ? n : 0);
  // structure.
  // the initial size should probably be exposed to the user
  struct cigar_string cigar = cigar_string_init(256);
  int qcig_length, q_beg, q_end;
  size_t ops_begin, ops_end;
  int rpos=0; // not actually used at the moment; needed as an argument to cigar_to_table
  struct cigar_parse_options cig_opt = init_cigar_parse_options();
  while( (*count) < n && read_next_entry(bam, b) >= 0){
    if( (f_flag & b->core.flag) != f_flag || (F_flag & b->core.flag) > 0 || b->core.qual < min_q )
      continue;
    if(ret_flag & (1 << S_QID))
      SET_STRING_ELT(VECTOR_ELT(ret_data, S_QID), *count, mkChar( bam_get_qname(b) ));
    if(ret_flag & (1 << S_FLAG))
      INTEGER(VECTOR_ELT(ret_data, S_FLAG))[*count] = b->core.flag;
    // add one to the ref identifier; this way it can be used with the return
    // value from the template lengths function.
    if(ret_flag & (1 << S_RNAME))
      INTEGER(VECTOR_ELT(ret_data, S_RNAME))[*count] = b->core.tid + 1;
    if(ret_flag & (1 << S_POS))
      INTEGER(VECTOR_ELT(ret_data, S_POS))[*count] = b->core.pos;
    if(ret_flag & (1 << S_MAPQ))
      INTEGER(VECTOR_ELT(ret_data, S_MAPQ))[*count] = b->core.qual;
    if(ret_flag & (1 << S_CIGAR)){
      cigar_string_read(&cigar, b);
      SET_STRING_ELT(VECTOR_ELT(ret_data, S_CIGAR), *count, mkChar( cigar.buffer ));
    }
    if(ret_flag & (1 << S_RNEXT))
      INTEGER(VECTOR_ELT(ret_data, S_RNEXT))[*count] = b->core.mtid + 1;
    if(ret_flag & (1 << S_PNEXT))
      INTEGER(VECTOR_ELT(ret_data, S_PNEXT))[*count] = b->core.mpos;
    if(ret_flag & (1 << S_TLEN))
      INTEGER(VECTOR_ELT(ret_data, S_TLEN))[*count] = b->core.isize;
    if(ret_flag & (1 << S_SEQ)){
      seq_buffer_size = bam_seq_p(b, &seq_buffer, seq_buffer_size);
      SET_STRING_ELT(VECTOR_ELT(ret_data, S_SEQ), *count, mkChar( seq_buffer ));
    }
    if(ret_flag & (1 << S_QUAL)){
      seq_buffer_size = bam_qual_p(b, &seq_buffer, seq_buffer_size);
      SET_STRING_ELT(VECTOR_ELT(ret_data, S_QUAL), *count, mkChar( seq_buffer ));
    }
    if(ret_flag & (1 << S_AUX)){
      // not currently making use of the return value. 
      bam_aux_string(b, &aux_str);
      SET_STRING_ELT(VECTOR_ELT(ret_data, S_AUX), *count, mkChar( aux_str.s ));
      aux_str.l = 0;
      aux_str.s[0] = 0;
    }
    if(ret_flag & (1 << S_CIG_TABLE)){
      cigar_to_table(b, 1+(*count), &rpos, &al_coord, &qcig_length, &q_beg, &q_end, &ops_begin, &ops_end, &cig_opt, 1);
      push_column(&query_info, (int[]){ b->core.l_qseq, qcig_length, q_beg, q_end,
					ops_begin, ops_end });
    }
    if(ret_flag & (1 << S_AUX_MM)){
      // the last two arguments provide the reference sequence and its length; we don't
      // have that here.
      parse_MM_string(b, 0, &mm_info, 1+(*count), &al_coord, ops_begin, ops_end, 0, 0);
    }
    (*count)++;
  }
  // if S_CIG_TABLE, make the appropriate data structures:
  if(ret_flag & (1 << S_CIG_TABLE)){
    SET_VECTOR_ELT( ret_data, S_CIG_TABLE, allocMatrix(INTSXP, al_coord.nrow, al_coord.col) );
    setAttrib(VECTOR_ELT(ret_data, S_CIG_TABLE), R_DimNamesSymbol, mk_rownames(cig_ops_rownames, CIG_OPS_RN));
    SET_VECTOR_ELT( ret_data, S_CIG_TABLE+1, allocMatrix(INTSXP, query_info.nrow, query_info.col) );
    setAttrib(VECTOR_ELT(ret_data, S_CIG_TABLE+1), R_DimNamesSymbol,
	      mk_rownames(q_info_rownames, Q_INFO_RN));
    memcpy( INTEGER(VECTOR_ELT(ret_data, S_CIG_TABLE)), al_coord.data, sizeof(int) * al_coord.nrow * al_coord.col );
    memcpy( INTEGER(VECTOR_ELT(ret_data, S_CIG_TABLE+1)), query_info.data, sizeof(int) * query_info.nrow * query_info.col );
  }
  if(ret_flag & (1 << S_AUX_MM)){
    SET_VECTOR_ELT( ret_data, S_AUX_MM, allocMatrix(INTSXP, mm_info.nrow, mm_info.col) );
    setAttrib(VECTOR_ELT(ret_data, S_AUX_MM), R_DimNamesSymbol, mk_rownames(mm_info_rownames, MM_INFO_RN));
    memcpy( INTEGER(VECTOR_ELT(ret_data, S_AUX_MM)), mm_info.data, sizeof(int) * mm_info.nrow * mm_info.col );
  }    
  free(seq_buffer);
  free(aux_str.s);
  free(al_coord.data);
  free(query_info.data);
  free(mm_info.data);
  cigar_string_free( &cigar );
  bam_destroy1(b);
  UNPROTECT(2);
  return(ret_data);
}

// Initially just return the basic information; that is no genotype data;
// I haven't yet quite worked out how the genotype data is encoded as it
// is not super clear from the data.. 
SEXP bcf_read_n(SEXP bcf_ptr_r, SEXP n_r, SEXP fmt_tags_r){
  struct bcf_ptrs *bcf = extract_bcf_ptr(bcf_ptr_r);
  if(!bcf || !bcf->vcf)
    error("External poiner or vcf/bcf file is NULL");
  if(TYPEOF(n_r) != INTSXP || length(n_r) != 1)
    error("n_r should be an integer vector of length 1");
  if(TYPEOF(fmt_tags_r) != STRSXP)
    error("fmt_tags_r should be a character vector");
  // bcf1_t, we have the following fields fields
  // pos, rlen:   64 bit signed integers (hts_pos_t, used to be 32 bits)
  // rid:         32 bit signed integer (int32_t)
  // qual:        float (32 bits)
  // n_info, n_allele:  unsigned 32 bit integers (but using 16 bits, so effectively shorts)
  //                    note that these use bitfields (which means that I don't think I can take
  //                    the address of them to store as a single integer.
  // n_fmt, n_sample:   unsigned 32 bit integers using 8 and 24 bits respectively.
  // kstring_t shared, indiv: Presumably a string holding the data;
  // use ks_c_str( kstring_t *s ) to make a null terminated const char
  //
  // the individual data is meant to be unpacked using bcf_unpack which will put it into a bcf_dec_t structure.
  // the chrom field is a signed 32 bit integer (int32_t)
  //
  // To specify a data structure for a set of FORMAT and INFO fields:
  // I need to have a mapping from [tag] -> [ret_data index].
  // This can be done via an intermediate array that has m values defaulting to -1
  // tag_index
  // tag_index[ tag_id ] -> ret_data index
  // the tag_id can be obtained using: bcf_hdr_id2int()
  // then the size of the tag_index is the maximum value of the tags requested.
  // But it may be simpler if I first define the formats for each tag, then
  // use, bcf_get_fmt_id (and bcf_get_info_id) for each record.
  // Work out the dimensions and types of sample specific data we will return:
  int *tag_ids = malloc(sizeof(int) * length(fmt_tags_r));
  int *tag_offsets = malloc(sizeof(int) * length(fmt_tags_r)); // return to indicate which tags identified
  int tag_n = 0;
  Rprintf("length of fmt_tags_r is %d\n", length(fmt_tags_r));
  for(int i=0; i < length(fmt_tags_r); ++i){
    int id = bcf_hdr_id2int(bcf->header, BCF_DT_ID, CHAR(STRING_ELT(fmt_tags_r, i)));
    if(id >= 0){
      tag_ids[tag_n] = id;
      tag_offsets[tag_n] = i;
      ++tag_n;
      Rprintf("%s: %ld\t%ld\t%d\t%d\t%d\n",
	      CHAR(STRING_ELT(fmt_tags_r, i)),
	      bcf_hdr_id2length( bcf->header, BCF_HL_FMT, id),
	      bcf_hdr_id2number( bcf->header, BCF_HL_FMT, id),
	      bcf_hdr_id2type( bcf->header, BCF_HL_FMT, id),
	      bcf_hdr_id2coltype( bcf->header, BCF_HL_FMT, id),
	      bcf_hdr_idinfo_exists( bcf->header, BCF_HL_FMT, id));
    }
  }

  // This would be better to define as a macro or const somewhere else.
  const char *const_field_names[12] = {"n", "rid", "rlen", "pos", "qual",
				       "n_info", "n_allele", "n_fmt", "n_sample",
				       "id", "als", "allleles"};
  int fields_n = 12 + tag_n;
  const char **field_names = malloc(sizeof(char*) * fields_n);
  for(int i=0; i < 12; ++i)
    field_names[i] = const_field_names[i];
  for(int i=0; i < tag_n; ++i)
    field_names[i + 12] = CHAR(STRING_ELT(fmt_tags_r, tag_offsets[i]));

  for(int i=0; i < 12 + tag_n; ++i){
    Rprintf("  %d : %s", i, field_names[i]);
  }
  Rprintf("\n");
  SEXP ret_data_names = PROTECT( mk_strsxp( field_names, fields_n));
					  
  /* SEXP ret_data_names = PROTECT( mk_strsxp( (const char*[]){"rid", "rlen", "pos", "qual", */
  /* 							      "n_info", "n_allele", "n_fmt", "n_sample", */
  /* 							      "shared", "indiv", "id", "als", "allleles", "n"}, 14)); */
  
  SEXP ret_data = PROTECT( allocVector(VECSXP, length(ret_data_names)) );
  setAttrib( ret_data, R_NamesSymbol, ret_data_names );
  // The first element is a count of the number of elements actually read.
  SET_VECTOR_ELT(ret_data, 0, allocVector(INTSXP, 1));
  int *count = INTEGER(VECTOR_ELT(ret_data, 0));
  *count = 0;
  int n=asInteger(n_r);
  if(n <= 0)
    error("you must request at least one entry");
  // allocate memory for the return data structure. Do this in a loop;
  const unsigned int BCF_ret_types[11] =
    {INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP,
     STRSXP, STRSXP, STRSXP};
  for(int i=0; i < 11; ++i)
    SET_VECTOR_ELT(ret_data, i+1, allocVector( BCF_ret_types[i], n ));
  
  bcf1_t *b = bcf_init();
  int which = BCF_UN_ALL; // unpack all information in order to work out how to use it.
  size_t buf_size = 256;
  char *ch_buf = malloc(sizeof(char) * buf_size);
  while((*count) < n && bcf_read(bcf->vcf, bcf->header, b) == 0){
    // unpack; this should fill the d (bcf_dec_t) data element
    bcf_unpack(b, which);
    // b->d has the following elements:
    // int m_fmt, m_info, m_id, m_als, m_allele, m_flt; ## these represent sizes of something
    // int n_flt; ## number of fields
    // int *flt; // FILTER keys in the dictionary ?
    // char *id, *als; // id and ref+alt block (0\separated). May be able to treat like const char
    // char **allele;  // allele[0] is the REF? Not sure how this should be treated.
    // bcf_info_t *info;
    // bcf_fmt_t *fmt;  // FORMAT and individual sample
    // and some more things. The actual allele information is probably in the fmt field.
    // extract some suitable information..
    INTEGER(VECTOR_ELT(ret_data, 1))[*count] = b->rid;
    REAL(VECTOR_ELT(ret_data, 2))[*count] = (double)b->rlen;
    REAL(VECTOR_ELT(ret_data, 3))[*count] = (double)b->pos;
    REAL(VECTOR_ELT(ret_data, 4))[*count] = (double)b->qual;
    INTEGER(VECTOR_ELT(ret_data, 5))[*count] = (int)b->n_info;
    INTEGER(VECTOR_ELT(ret_data, 6))[*count] = (int)b->n_allele;
    INTEGER(VECTOR_ELT(ret_data, 7))[*count] = (int)b->n_fmt;
    INTEGER(VECTOR_ELT(ret_data, 8))[*count] = (int)b->n_sample;
    SET_STRING_ELT(VECTOR_ELT(ret_data, 9), *count, mkChar(b->d.id));
    SET_STRING_ELT(VECTOR_ELT(ret_data, 10), *count, mkChar(b->d.als));
    // set only the reference allele to start with. Just in order to test how to deal with the
    // data
    size_t buf_offset = 0;
    // the following code combines a set of 0 terminated chars. There should be a function
    // for this; either in this code or in string.h or something.
    for(int i=0; i < b->n_allele; ++ i){
      int len = strlen( b->d.allele[i] );
      if( 1 + buf_offset + len > buf_size ){
	buf_size = 2 * (1 + buf_offset + len);
	ch_buf = realloc(ch_buf, sizeof(char) * buf_size );
	// in theory, ch_buf could be 0 at this point, if it is, we should crash shortly.
      }
      memcpy( ch_buf + buf_offset, b->d.allele[i], len );
      ch_buf[ buf_offset + len ] = ',';
      buf_offset += (len + 1);
    }
    ch_buf[ buf_offset - 1] = 0; 
    SET_STRING_ELT(VECTOR_ELT(ret_data, 11), *count, mkChar(ch_buf));
    //    SET_STRING_ELT(VECTOR_ELT(ret_data, 12), *count, mkChar(b->d.allele[0]));
    //// TEMPORARY BLOCK FOR TESTING //////
    const char *fm_id[3] = {"GT", "GQ", "AD"};
    for(int i=0; i < 3; ++i){
      bcf_fmt_t *fmt = bcf_get_fmt( bcf->header, b, fm_id[i] );
      if(fmt != 0){
	Rprintf("%s : %d  %d %d\n", fm_id[i], fmt->n, fmt->size, fmt->type );
      }else{
	Rprintf("%s : no entry found \n", fm_id[i]);
      }
    }
    for(int i=0; i < tag_n; ++i){
      Rprintf("%s : \n", CHAR(STRING_ELT(fmt_tags_r, tag_offsets[i])));
      // the call to bcf_get_format_values crashes here. Maybe we should just get the data directly using the fmt information.
      //      int nn = bcf_get_format_values( bcf->header, b, CHAR(STRING_ELT(fmt_tags_r, tag_offsets[i])), (void*)dst, &ndst, BCF_HT_INT);
      bcf_fmt_t *fmt = bcf_get_fmt_id(b, tag_ids[i]);
      // we have fmt.n, fmt.size, fmt.type; n is the number of values per sample. size the number of bytes per sample and type the data type.
      Rprintf("%d  %d  %d  %d  %d\n", fmt->id, fmt->n, fmt->size, fmt->type, fmt->p_len);
      int ns = fmt->p_len / fmt->size;
      // OK; this is a major pain. But it does work... 
      for(int j=0; j < ns && j < 5; ++j){
	for(int k=0; k < fmt->n; ++k){
	  int off = j * fmt->size + k * (fmt->size / fmt->n);
	  Rprintf("%d %d %d:  ", j, k, off);
	  switch(fmt->type){
	  case BCF_BT_INT8:
	    Rprintf(" %d ", (int)*(fmt->p + off));
	    break;
	  case BCF_BT_INT16:
	    Rprintf(" %d ", (int)*((short*)fmt->p + off));
	    break;
	  case BCF_BT_INT32:
	    Rprintf(" %d ", *((int*)fmt->p + off));
	    break;
	  case BCF_BT_FLOAT:
	    Rprintf(" %f ", *((float*)fmt->p + off));
	    break;
	  default:
	    Rprintf("don't know how to do that\n");
	  }
	  Rprintf("\n");
	}
      }
    }
    
    ////////////////////////////////////////////
    //    Rprintf("%d\t%d\n", b->n_fmt, b->d.m_fmt);
    // d.m_fmt seems to be the maximum number of fmts;
    // m is the size, n is the number.. we should be using n everywhere we have it.. 
    /* for(int i=0; i < b->n_fmt; ++i) */
    /*   Rprintf("\t%d: %s", b->d.fmt[i].id, bcf->header->id[BCF_DT_ID][ b->d.fmt[i].id ].key); */
    /* Rprintf("\n"); */
    ++(*count);
  }
  free(ch_buf);
  free(tag_ids);
  free(tag_offsets);
  free(field_names);
  bcf_destroy1(b);
  UNPROTECT(2);
  return(ret_data);
}

// what to return depends on the bitwise ret_flag_r
// 01:   return individual bit counts (vector of size 12)
// 10:   return a vector of all 4096 different flag values
// Somewhat counterintuitively the latter is likely to be
// Regardless, two vectors are returned in any case.
SEXP sam_flag_stats(SEXP bam_ptr_r, SEXP ret_flag_r){
  struct bam_ptrs *bam = extract_bam_ptr(bam_ptr_r);
  if(!bam || !bam->sam){
    // we could be clever here and try to recreate from the protect information
    // but we should think carefully about this
    error("External pointer or sam file is NULL");
  }
  if(TYPEOF(ret_flag_r) != INTSXP || length(ret_flag_r) != 1)
    error("ret_flag_r should be an integer vector of length 1");
  unsigned int ret_flag = (unsigned int)asInteger(ret_flag_r);
  if(!(ret_flag & 3))
    error("Neither bit 1 or 2 set in ret_flag:");
  // The sam defines 12 different bits; we can return a vector
  // of the individual bits and a count of all of the complete combinations
  unsigned int sam_bit_n = 12;
  unsigned int mask = (1 << sam_bit_n) - 1;
  int *bit_count = 0;
  int *flag_count = 0;
  SEXP ret_data = PROTECT(allocVector(VECSXP, 2));
  if(ret_flag & 1){
    SET_VECTOR_ELT( ret_data, 0, allocVector(INTSXP, sam_bit_n));
    bit_count = INTEGER(VECTOR_ELT( ret_data, 0));
    memset(bit_count, 0, sizeof(int) * sam_bit_n);
  }
  if(ret_flag & 2){
    SET_VECTOR_ELT( ret_data, 1, allocVector(INTSXP, 1 << sam_bit_n));
    flag_count = INTEGER(VECTOR_ELT( ret_data, 1));
    memset(flag_count, 0, sizeof(int) * (1 << sam_bit_n));
  }
  // and then go through and count the bits:
  bam1_t *b = bam_init1();
  while( sam_read1(bam->sam, bam->header, b) >= 0 ){
    unsigned int flag = b->core.flag & mask;
    if(ret_flag & 1){
      for(unsigned int j=0; j < sam_bit_n; ++j)
	bit_count[j] += (flag & (1 << j)) ? 1 : 0;
    }
    if(ret_flag & 2)
      flag_count[flag]++;
  }
  bam_destroy1(b);
  UNPROTECT(1);
  return(ret_data);
}

// simply return BAM_CIGAR_STR as a STRSXP object
SEXP bam_cigar_str(){
  SEXP bcs = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT( bcs, 0, mkChar(BAM_CIGAR_STR) );
  UNPROTECT(1);
  return(bcs);
}

// simply return BAM_CIGAR_STR as a STRSXP object
SEXP nuc_table(){
  SEXP nuc_r = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT( nuc_r, 0, mkChar(nuc_encoding) );
  UNPROTECT(1);
  return(nuc_r);
}

// flags..
SEXP bam_flag(SEXP flags_r){
  if(TYPEOF(flags_r) != INTSXP || length(flags_r) < 1)
    error("bam_flag flags_r must be an integer vector of positive length");
  size_t l = length(flags_r);
  SEXP ret_data = PROTECT(allocVector(STRSXP, l));
  int *flags = INTEGER(flags_r);
  for(int i=0; i < l; ++i){
    char *flag_string = bam_flag2str(flags[i]);
    SET_STRING_ELT(ret_data, i, mkChar(flag_string));
    free(flag_string);
  }
  UNPROTECT(1);
  return(ret_data);
}
  


// Return the lengths of the sequences
// as a named vector
// CHANGE: have a function to validate the bam_ptr_r as we use that more than
// once.
SEXP target_lengths(SEXP bam_ptr_r){
  struct bam_ptrs *bam = extract_bam_ptr(bam_ptr_r);
  if(!bam){
    // we could be clever here and try to recreate from the protect information
    // but we should think carefully about this
    error("External pointer is NULL");
  }
  SEXP t_lengths_r = PROTECT(allocVector(INTSXP, bam->header->n_targets));
  SEXP t_names_r = PROTECT(allocVector(STRSXP, bam->header->n_targets));
  int *t_lengths = INTEGER(t_lengths_r);
  for(int i=0; i < bam->header->n_targets; ++i){
    t_lengths[i] = bam->header->target_len[i];
    SET_STRING_ELT( t_names_r, i, mkChar(bam->header->target_name[i]) );
  }
  setAttrib( t_lengths_r, R_NamesSymbol, t_names_r );
  UNPROTECT(2);
  return(t_lengths_r);
}

// This should only be called from extract_aux_tags
// As the terms are very specific; 
void set_aux_default_values(unsigned int i, char *aux_types, int n, SEXP ret_data, uint64_t tags_set){
  for(int k=0; k < n; ++k){
    if(((tags_set >> k) & 1) == 1)
      continue;
    switch(aux_types[k]){
    case 'A':
    case 'C':
    case 'S':
    case 'I':
      INTEGER(VECTOR_ELT(ret_data, k))[i] = R_NaInt;
      break;
    case 'F':
      REAL(VECTOR_ELT(ret_data, k))[i] = R_NaN;
      break;
    default:
      SET_STRING_ELT( VECTOR_ELT(ret_data, k), i, mkChar(""));
    }
  }
}

// This function extracts the values of specific auxiliar fields
// from aux strings. This could probably be done more efficiently
// by using the bam representation and calling the accessor functions
// defined in htslib; and included in sam.read.n. However, that would
// require careful reading of the specific libraries again, and it
// could not be used with arbitrary aux strings obtained from elsewhere.
// This will return a list with one member for each tag type specified.
// Each list will be the same length as aux_r, the number of auxiliar strings.
SEXP extract_aux_tags(SEXP aux_r, SEXP aux_tags_r, SEXP aux_types_r){
  if(TYPEOF(aux_r) != STRSXP || length(aux_r) < 1)
    error("aux_r must be a string with length of 1 or longer");
  if(TYPEOF(aux_tags_r) != STRSXP || length(aux_tags_r) < 1)
    error("aux_tags_r must be a string with length of 1 or longer");
  if(TYPEOF(aux_types_r) != STRSXP || length(aux_types_r) < 1 || length(aux_types_r) > 64)
    error("aux_types_r must be a string with length between 1 and 64 (inclusive)");
  if(length(aux_types_r) != length(aux_tags_r))
    error("The length of aux_types_r and aux_tags_r must be the same");
  // tag values should be two characters;
  // check that this is correct before anything doing anything else:
  for(int i=0; i < length(aux_tags_r); ++i){
    if(length(STRING_ELT(aux_tags_r, i)) != 2){
      Rprintf("Illegal tag: %s\n", CHAR(STRING_ELT(aux_tags_r, i)));
      error("Illegal tag specified");
    }
  }
  // the types are represented by a single character which
  // must be one of: A, c, C, s, S, i, I, f, Z, H
  // where: A single char, c int8, C uint 8, s int16, S uint16, i int32, I uint32, f float (32 bit)
  // Z: null terminated char array, H: null terminated char array (of Hex representations)
  // We can distinguish the upper case of all of these by their last 5 bits
  // which give values of, 1, 3, 19, 9, 6, 26, and 8 respectively.
  // That means we can check if a type is defined by AND operations on a 32 bit unsigned
  // integer
  // note that toupper is equivalent to & 0xDF
  uint8_t up_mask = 0xDF;
  // This is with bits 1, 3, 19, 9, 6, 26 and 8 set.
  uint32_t type_check = 0x020401A5;
  // to access the upper three bits of a char:
  uint8_t left_mask = 0xE0;
  uint8_t left_check = 0x40;
  // this will give the 5 lower bits
  // that can be used as bit offset that can be checked with the type_check.
  uint8_t right_mask = 0x1F;
  // first go through and check that the types specified are correct. Then go through and
  // check if the tags are correct (each should be a two letter code).
  // use an array 
  char *aux_types = malloc(length(aux_types_r));
  for(int i=0; i < length(aux_types_r); ++i){
    SEXP type_r = STRING_ELT(aux_types_r, i);
    if(length(type_r) != 1){
      free(aux_types);
      error("All type specifiers must be single characters");
    }
    aux_types[i] = CHAR(type_r)[0] & up_mask; // sets to upper as we don't distinguish
    // check if it is allowed;
    if( !(((1 << ((aux_types[i] & right_mask) - 1)) & type_check) && ((left_mask & aux_types[i]) == left_check) ) ){
      Rprintf("illegal aux_type: %c\n", aux_types[i]);
      free(aux_types);
      error("Error illegal aux_type");
    }
  }
  // set up some data structures:
  SEXP ret_data = PROTECT(allocVector(VECSXP, length(aux_types_r)));
  for(int i=0; i < length(ret_data); ++i){
    switch(aux_types[i]){
    case 'A':
    case 'C':
    case 'S':
    case 'I':
      SET_VECTOR_ELT(ret_data, i, allocVector(INTSXP, length(aux_r)));
      break;
    case 'F':
      SET_VECTOR_ELT(ret_data, i, allocVector(REALSXP, length(aux_r)));
      break;
    default: // this is a pain
      SET_VECTOR_ELT(ret_data, i, allocVector(STRSXP, length(aux_r)));
    }
  }
  // then go through and allocate the appropriate type of memory; int, double or string
  // if 
  // Then go through each aux_string; split on "\t", check if the tag is defined
  for(int i=0; i < length(aux_r); ++i){
    const char *aux = CHAR(STRING_ELT(aux_r, i));
    // This is likely slow since it will make lots of new words.
    // And will hence allocate lots of new strings. We can make it faster
    // quite easily by simply defining the offsets and lengths of all terms
    // instead. But first make it work.
    unsigned int n = 0;
    char **terms = split_string(aux, '\t', &n);
    int *int_values;
    double *double_values;
    // WARNING: the use of tags_set means that we are limiting ourselves to
    // a maximum of 64 tags. I need to check that somewhere above;
    // if a value for a tag is found the corresponding bit will be set to 1
    uint64_t tags_set = 0;
    for(unsigned int j=0; j < n; ++j){
      // check and then free the word..
      // the term should be a minimum of length 6
      // with : at positions 3 and 5 (counting from 1)
      size_t tl = strlen(terms[j]);
      if(tl >= 6 && terms[j][2] == ':' && terms[j][4] == ':'){
	// this is slow as well
	for(int k=0; k < length(aux_tags_r); ++k){
	  if(strncmp(CHAR(STRING_ELT(aux_tags_r, k)), terms[j], 2) == 0){
	    // do something useful depending on if the type
	    // fits.
	    if(aux_types[k] == (terms[j][3] & up_mask)){
	      tags_set |= (1 << k);
	      // Everything fits. Now we need to do something different depending on what the
	      // type is.
	      switch(aux_types[k]){
	      case 'A':
		int_values = INTEGER(VECTOR_ELT(ret_data, k));
		int_values[i] = (int)terms[j][5];
		break;
	      case 'C':
	      case 'S':
	      case 'I':
		int_values = INTEGER(VECTOR_ELT(ret_data, k));
		int_values[i] = atoi( terms[j] + 5 );
		break;
	      case 'F':
		double_values = REAL(VECTOR_ELT(ret_data, k));
		double_values[i] = atof( terms[j] + 5 );
		break;
	      default:
		SET_STRING_ELT( VECTOR_ELT(ret_data, k), i, mkChar(terms[j] + 5));
	      }
	    }
	    break;
	  }
	}
      }
      // need a function to set default values..
      set_aux_default_values(i, aux_types, length(aux_tags_r), ret_data, tags_set);
      free( terms[j] );
    }
    free(terms);
  }
  free(aux_types);
  // case then go through and check the values.
  UNPROTECT(1);
  return(ret_data);
}

// Make use of query_to_ref to map from positions in unaligned query
// sequences to to reference positions.
// Arguments:
// q_pos: a table containing the following columns:
//   1. query positions; [sorted or not?]
//   2. first row in an ops table for the alignment
//   3. last row in an ops table for the alignment
//   4. an integer value specifying direction (0 = reverse)
//   5. The query length
// ops: an ops table containing usual rows of the usual number of rows:
//   1. al.i
//   2. op
//   3. type
//   4. r0
//   5. q0
//   6. r1
//   7. q1
//   8. op.l
// Not all these values will be used, but query_to_ref expects an i_matrix of these
// dimensions. And it's easier to make one in R than to do it here.
// Specifically one can simply return a table with only match operations (these need
// to be 0)
//
// Returns: For each query position
// A matrix containing the following columns:
// 1. query position (a copy)
// 2. reference position
// 3. the index of the cigar operation within which the match was found.

SEXP qpos_to_ref_pos(SEXP qpos_r, SEXP ops_r){
  if(TYPEOF(qpos_r) != INTSXP || TYPEOF(ops_r) != INTSXP){
    //    warning("qpos_to_ref_pos: both arguments should be integer matrices");
    return(R_NilValue);
  }
  SEXP qpos_dim_r = getAttrib(qpos_r, R_DimSymbol);
  SEXP ops_dim_r =  getAttrib(ops_r, R_DimSymbol);
  if(length(qpos_dim_r) != 2 || length(ops_dim_r) != 2){
    //    warning("qpos_to_ref_pos: both arguments should be matrices");
    return(R_NilValue);
  }
  int *qpos_dim = INTEGER(qpos_dim_r);
  int *ops_dim = INTEGER(ops_dim_r);
  if(qpos_dim[0] < 1 || qpos_dim[1] != 5){
    //    warning("The query positions table must have at least one row and exactly five columns");
    return(R_NilValue);
  }
  if(ops_dim[0] != CIG_OPS_RN || ops_dim[1] < 1){
    //    warning("The ops table should have exactly 8 rows and at least one column");
    return(R_NilValue);
  }
  struct i_matrix ops = init_i_matrix(ops_dim[0], ops_dim[1]);
  memcpy(ops.data, INTEGER(ops_r), sizeof(int) * ops_dim[0] * ops_dim[1]);
  ops.col=ops_dim[1];
  int *qpos_table = INTEGER(qpos_r);
  int *qpos = qpos_table;
  int *ops_begin = qpos_table + qpos_dim[0];
  int *ops_end = ops_begin + qpos_dim[0];
  int *is_fwd = ops_end + qpos_dim[0];
  int *qlen = is_fwd + qpos_dim[0];
  int r_pos = -1;
  SEXP ret_data_r = PROTECT( allocMatrix(INTSXP, qpos_dim[0], 3) );
  int *ret_data = INTEGER(ret_data_r);
  int *ret_qpos = ret_data;
  int *ret_rpos = ret_data + qpos_dim[0];
  int *ret_op_i = ret_rpos + qpos_dim[0];
  int in_clipped = 2; // try to assign to the adjacent region if position is not in an M operation
  for(int i=0; i < qpos_dim[0]; ++i){
    if(ops_begin[i] < 0 || ops_begin[i] > ops_end[i] || ops_end[i] > ops_dim[1]){
      //      warning("Inappropriate ops_begin or ops_end values: %d  %d skipping", ops_begin[i], ops_end[i]);
      continue;
    }
    // ops_begin[i]-1, because the table has been passed from R.
    // in R, the end is usually the last index, rather than index + 1
    // Hence we do not need to subtract from ops_end[i]
    size_t col_i = is_fwd[i] ? (size_t)ops_begin[i]-1 : (size_t)ops_end[i]-1;
    int r = query_to_ref(qpos[i], &r_pos, &ops, ops_begin[i]-1, ops_end[i],
			 &col_i, is_fwd[i], qlen[i], in_clipped);
    ret_qpos[i] = r;
    ret_rpos[i] = r_pos;
    ret_op_i[i] = (int)col_i + 1;
  }
  free( ops.data );
  UNPROTECT(1);
  return( ret_data_r );
}

/////// This function got too complicated; see option below for alternative strategy ///
// merge_transcripts
// given an ops table with columns:
//   al.i, op, type, r0, q0, r1, q1, op.l
// a vector of flags
// a maximum 3' distance
// a vector of [required flags, banned flags ]
//    Note that one and only one of these must have 16 set
//    i.e. reg_f XOR ban_f must equal 16.
// Will merge alignments assumed to be long transcripts following
// some simple rules:
// 1. The 3' end of the alignments must be within the max distance
//    (this could be considered optional; the main issue isn't what
//    is correct biologically, but what is easy to encode.
// 2. All upstream intron positions must be identical.
//
// This doesn't really feel like it belongs here, but since it is making
// use of the ops table structure it sort of makes sense to use it here.
/* SEXP merge_transcripts(SEXP ops_r, SEXP flags_r, */
/* 		       SEXP max_distance_r, SEXP flag_filter_r){ */
/*   // In order to simplify the logic, this function will make use of a set */
/*   // of locally defined structs; these will contain only integers; */

/*   typedef struct transcripts { */
/*     // p_start and p_end are chromosomal coordinates giving the range */
/*     // p_start < p_end */
/*     // row_start and row_end are indices in the transcript_rows table */
/*     int p_start; int p_end; int row_start; int row_end; */
/*   } transcripts; */
/*   typedef struct transcript_rows { */
/*     // tr_i: the index in a transcripts table */
/*     // ops_row: a row in an cigar_ops table */
/*     // count: the number of times a the operation has been encountered */
/*     // next_step: the distance to the next entry in this table */
/*     //            0 if no further entries for this transcript */
/*     int tr_i; int ops_row; int count; int next_step; */
/*   } transcript_rows; */
/*   typedef struct transcript_alignment { */
/*     // tr_i: the index of an entry in transcripts */
/*     // al_i: the identifier of an alignment. Taken directly from the ops */
/*     //       table. */
/*     int tr_i; int al_i; */
/*   } transcript_alignment; */
  
/*   if(TYPEOF(ops_r) != INTSXP || length(ops_r) < 8) */
/*     error("ops_r should be an integer matrix with 8 columns"); */
/*   if(TYPEOF(flags_r) != INTSXP || TYPEOF(max_distance_r) != INTSXP) */
/*     error("flags_r AND max_distance_r should be integer vectors"); */
/*   if(TYPEOF(flag_filter_r) != INTSXP || length(flag_filter_r) != 2) */
/*     error("flag_filter_r should be an integer vector of length 2"); */
/*   SEXP op_dims_r = getAttrib(ops_r, R_DimSymbol); */
/*   if(length(op_dims_r) != 2) */
/*     error("ops_r should be a matrix"); */
/*   int *dims = INTEGER(op_dims_r); */
/*   // As of writing we must have at least 8 columns (dims[1]) */
/*   int ops_ncol = 8; */
/*   if(dims[1] != ops_ncol) */
/*     error("The ops_r table should have 8 columns. Specified table has: %d", dims[1]); */
/*   int ops_nrow = dims[0]; */
/*   if(ops_nrow < 1) */
/*     error("The ops_r table must have at least one row"); */
/*   int *ops = INTEGER(ops_r); */
/*   int *ops_ali = ops; */
/*   int *ops_op = ops + ops_nrow; */
/*   int *ops_r0 = ops + (3 * ops_nrow); */
/*   int *ops_r1 = ops + (5 * ops_nrow); */

/*   // check that the flags vector is of the correct length: */
/*   // This assumes that the ops is sorted by alignment id with the largest */
/*   // last: */
/*   // and -1 because the table is from R where counting is from 1. */
/*   int flags_n = length(flags_r); */
/*   if(flags_n != (ops_ali[ops_nrow-1] - 1)) */
/*     error("There must be a flag for every entry"); */
/*   int *flags = INTEGER(flags_r); */

/*   if(length(max_distance_r) != 1) */
/*     error("max_distance_r should have a length of 1"); */
/*   int max_distance = INTEGER(max_distance_r)[0]; */
  
/*   // The flag_filter must specify the orientation of reads; forward or reverse */
/*   // This implies that flag 16 (0x10) must either be required or banned */
/*   // It is also clear that no flags should both be required and banned */
/*   int *flag_filter = INTEGER(flag_filter_r); */
/*   // As I understand it, the following operation is technically undefined */
/*   // for signed integers; but R doesn't have unsigned ints... */
/*   if(flag_filter[0] & flag_filter[1] > 0 || */
/*      (0x10 & (flag_filter[0] ^ flag_filter[1])) != 0x10) */
/*     error("Improper flag filter set; either common flags or direction not specified"); */
  
/*   // To store the transcript information we make use of four tables: */
/*   // tr: p_start, p_end, row_start, row_end */
/*   //    p_start: the genomic start position (p_start < p_end) */
/*   //    p_end:   the genomic end position (p_end > p_start) */
/*   //    row_start: first row in tr_rows */
/*   //    row_end:   last row in tr_rows */
/*   //  */
/*   // tr_rows: tr_i, row, count, step_next */
/*   // The rows in the ops table that define a transcript */
/*   //    tr_i: 0 based index of the transcript. These refer to the row in the tr table. */
/*   //    count: the number of times that the op has been seen in compatible */
/*   //           transcripts */
/*   //    step_next: the number of rows in tr_rows before the appearance of the next */
/*   //               op for the given transcript. This is necessary because compatible */
/*   //               transcripts may not appear consecutively in the alignments. */
/*   //  note; an entry is made only the first time a unique operation is encountered; */
/*   //        after that, the count is incremented. */
/*   //  */
/*   // tr_al: tr_i, al_i, */
/*   // where: tr_i is a transcript identifier; */
/*   //        al_i are the alignment identifiers from the ops table */
/*   // AND we need to keep track of the smallest row in the tr_pos table */
/*   // that is still active. */

/*   // note that: */
/*   // tr_vec.length will give the number of transcrips that have been defined */
/*   struct tr_vec = vector_init( OPS_INIT_SIZE, sizeof(transcripts) ); */
/*   struct tr_rows_vec = vector_init( OPS_INIT_SIZE, sizeof(transcript_rows) ); */
/*   struct tr_al_vec = vector_init( OPS_INIT_SIZE, sizeof(transcript_alignment) ); */

/*   // We will also keep track of matching rows in tr_rows_vec */
/*   struct vectori tr_matching_rows = vectori_init( OPS_INIT_SIZE ); */
/*   struct vectori al_matching_rows = vectori_init( OPS_INIT_SIZE ); */
/*   // the oldest transcript that can still overlap with the last */
/*   // observed alignment. */
/*   int oldest_active_tr = 0; */
/*   int current_al_i = ops_ali[0]; */
/*   //  int last_al_i = ops_ali[0]; */
/*   int ops_start_row = 0; */
/*   int ops_end_row = 0; */
/*   for(int i=0; i < ops_nrow; ++i){ */
/*     last_al_i = current_al_i = ops_ali[i]; */
/*     // check flags; */
/*     // and skip until we get to an ok flag and a match operation */
/*     while( i < ops_nrow && */
/* 	   flag_filter[0] & flags[ ops_ali[i]-1 ] == 0 && */
/* 	   flag_filter[1] & flags[ ops_ali[i]-1 ] > 0 && */
/* 	   ops_ops[i] != CIG_M){ */
/*       ++i; */
/*     } */
/*     if(i == ops_nrow) */
/*       break; */
/*     current_al_i = ops_ali[i]-1; */
/*     ops_start_row = i; */
/*     int al_start = ops_r0[i]; */
/*     // then find the end of the current alignment: */
/*     int al_end = -1; */
/*     while(i < ops_nrow && ops_ali[i] == current_al_i){ */
/*       if(ops_op[i] == CIG_M){ */
/* 	al_end = ops_r1[i]; */
/* 	ops_end_row = i; */
/*       } */
/*       ++i; */
/*     } */
/*     if(al_end == -1){ */
/*       Rprintf("failed to find end for alignment: %d\n, skipping\n", current_al_i + 1); */
/*       continue; */
/*     } */
/*     // Then go through all active alignments and check compatibility; */
/*     int compat_tr_n = 0; */
/*     // alignment is forward if 0x10 is banned...  */
/*     int fwd = flag_filter[1] & 0x10 > 0; */
/*     transcript *tr = 0; */
/*     for(int j=oldest_active_tr; j < tr_vec.length; ++j){ */
/*       tr = (transcript *tr)vector_at(&tr, j); */
/*       if(tr->p_end < ops_r0[al_start]){ */
/* 	oldest_active_tr = j + 1; */
/* 	continue; */
/*       } */
/*       // check for compatability; */
/*       int k = tr->row_start; */
/*       int l = ops_start_row;  // must stay within ops_start_row -> ops_end_row */
/*       int compatible = 1; */
/*       // clear the record of matching rows; */
/*       tr_matching_rows.n = 0; */
/*       al_matching_rows.n = 0; */
/*       // we can skip directly if fwd is FALSE and the start points are too divergent */
/*       if(!fwd && abs( al_start - tr->p_start ) > max_distance) */
/* 	continue; */
/*       while(k < tr_rows_vec.length){ */
/* 	transcript_rows *tr_rows = (transcript_rows*)vector_at(&tr_rows_vec, k); */
/* 	int increment = tr_rows.next_step; */
/* 	// increment l until we hit an intron; */
/* 	while( l < ops_end_row && ops_op[l] != CIG_N ) */
/* 	  ++l; */
/* 	// if the operation in tr is not an intron, go to next */
/* 	if( ops_op[ tr_rows->ops_row ] != CIG_N ) */
/* 	  goto next_op; */

/* 	int introns_match = (ops_r0[ tr_rows->ops_row ] == ops_r0[l] && ops_r1[ tr_rows->ops_row ] == ops_r1[l]); */
/* 	// if !fwd the intron operations must match; */
/* 	if(!fwd && !introns_match){ */
/* 	  compatible = 0; */
/* 	  break; */
/* 	} */
/* 	// if fwd, the intron operations must match once the tr intron is caugth up with */
/* 	while(l < ops_end_row && ops_op[l] != CIG_N && ops_r0[l] < ops_r0[ tr_rows->ops_row ]) */
/* 	  ++l; */
/* 	introns_match = (ops_r0[ tr_rows->ops_row ] == ops_r0[l] && ops_r1[ tr_rows->ops_row ] == ops_r1[l]); */
/* 	if(!introns_match){ */
/* 	  compatible = 0; */
/* 	  break; */
/* 	} */
/* 	vectori_push( tr_matching_rows, k ); */
/* 	vectori_push( al_matching_rows, l ); */
/* 	if(increment == 0) // last entry; should have pointed to a match */
/* 	  break; */
/*       next_op: */
/* 	k += increment; */
/*       } */
/*       // if fwd, then the last tr operation must m */
/*       if(!compatible || (fwd && abs( al_end - tr->p_end ) > max_distance)) */
/* 	continue; */
/*       // If we get here, then compatible should be true; */
/*       for(k=0; mi < tr_matching_rows.n; ++mi){ */
/* 	tr_row = (transcript_rows*)vector_at(&tr_rows_vec, tr_matching_rows.data[k]); */
/* 	tr_row->count++; */
/*       } */
/*       // if fwd is true, then either: */
/*       // increment the counter for the first (match) operation or: */
/*       // add an additional match operation at the end. This will need to have */
/*       if(fwd){ */
/* 	int k = tr->row_start; */
/* 	tr_rows = (transcript_rows*)vector_at(&tr_rows_vec, k); */
/* 	transcripts_rows *prev_tr_row = tr_rows; */
/* 	int increment; */
/* 	while(1){ */
/* 	  increment = tr_rows->next_step; */
/* 	  if(increment == 0) */
/* 	    break; */
/* 	  tr_rows = (transcript_rows*)vector_at(&tr_rows_vec, k + increment); */
/* 	  if(ops_op[ tr_rows->ops_row ] != CIG_M) */
/* 	    break; */
/* 	  k += increment; */
/* 	} */
/* 	if( ops_r0[ prev_tr_row->ops_row ] == tr->p_start ){ */
/* 	  prev_tr_row->count++; */
/* 	}else{ */
/* 	  // do complicated merging procedure... */
/* 	  int new_step = k - tr_rows_vec.length + increment; */
/* 	  prev_tr_row->next_step = tr_rows_vec.length - k; */
/* 	  vector_push( &tr_rows_vec, (void*){ prev_tr_row->tr-i, prev_tr_row->ops_row, 1, new_step }); */
/* 	} */
/*       } */
/*       // if not forward we need to add all of the additional rows from l to end;  */
/*     } */
/*   } */

/* } */


// merge_transcripts
// given an ops table with columns:
//   al.i, op, type, r0, q0, r1, q1, op.l
// a vector of flags
// a maximum 3' distance
// a vector of [required flags, banned flags ]
//    Note that one and only one of these must have 16 set
//    i.e. reg_f XOR ban_f must equal 16.
// Will merge alignments assumed to be long transcripts following
// some simple rules:
// 1. The 3' end of the alignments must be within the max distance
//    (this could be considered optional; the main issue isn't what
//    is correct biologically, but what is easy to encode.
// 2. All upstream intron positions must be identical.


// This should be refactored to use code in a separate compilation unit
// This one is getting too complex.
/* SEXP merge_transcripts(SEXP ops_r, SEXP flags_r, */
/* 		       SEXP max_distance_r, SEXP flag_filter_r){ */
/*   // In order to simplify the logic, this function will make use of a set */
/*   // of locally defined structs; these will contain only integers; */
/*   if(TYPEOF(ops_r) != INTSXP || length(ops_r) < 8) */
/*     error("ops_r should be an integer matrix with 8 columns"); */
/*   if(TYPEOF(flags_r) != INTSXP || TYPEOF(max_distance_r) != INTSXP) */
/*     error("flags_r AND max_distance_r should be integer vectors"); */
/*   if(TYPEOF(flag_filter_r) != INTSXP || length(flag_filter_r) != 2) */
/*     error("flag_filter_r should be an integer vector of length 2"); */
/*   SEXP op_dims_r = getAttrib(ops_r, R_DimSymbol); */
/*   if(length(op_dims_r) != 2) */
/*     error("ops_r should be a matrix"); */
/*   int *dims = INTEGER(op_dims_r); */
/*   // As of writing we must have at least 8 columns (dims[1]) */
/*   int ops_ncol = 8; */
/*   if(dims[1] != ops_ncol) */
/*     error("The ops_r table should have 8 columns. Specified table has: %d", dims[1]); */
/*   int ops_nrow = dims[0]; */
/*   if(ops_nrow < 1) */
/*     error("The ops_r table must have at least one row"); */
/*   int *ops = INTEGER(ops_r); */
/*   int *ops_ali = ops; */
/*   int *ops_op = ops + ops_nrow; */
/*   int *ops_r0 = ops + (3 * ops_nrow); */
/*   int *ops_r1 = ops + (5 * ops_nrow); */
  
/*   // check that the flags vector is of the correct length: */
/*   // This assumes that the ops is sorted by alignment id with the largest */
/*   // last: */
/*   // and -1 because the table is from R where counting is from 1. */
/*   int flags_n = length(flags_r); */
/*   if(flags_n != (ops_ali[ops_nrow-1] - 1)) */
/*     error("There must be a flag for every entry"); */
/*   int *flags = INTEGER(flags_r); */
  
/*   if(length(max_distance_r) != 1) */
/*     error("max_distance_r should have a length of 1"); */
/*   int max_distance = INTEGER(max_distance_r)[0]; */
  
/*   // The flag_filter must specify the orientation of reads; forward or reverse */
/*   // This implies that flag 16 (0x10) must either be required or banned */
/*   // It is also clear that no flags should both be required and banned */
/*   int *flag_filter = INTEGER(flag_filter_r); */
/*   // As I understand it, the following operation is technically undefined */
/*   // for signed integers; but R doesn't have unsigned ints... */
/*   if(flag_filter[0] & flag_filter[1] > 0 || */
/*      (0x10 & (flag_filter[0] ^ flag_filter[1])) != 0x10) */
/*      error("Improper flag filter set; either common flags or direction not specified"); */
/* } */

static const R_CallMethodDef callMethods[] = {
  {"load_bam", (DL_FUNC)&load_bam, 2},
  {"set_iterator", (DL_FUNC)&set_iterator, 3},
  {"clear_iterator", (DL_FUNC)&clear_iterator, 1},
  {"build_query_index", (DL_FUNC)&build_query_index, 4},
  {"query_positions", (DL_FUNC)&query_positions, 2},
  {"alignments_region", (DL_FUNC)&alignments_region, 9},
  {"alignments_region_mt", (DL_FUNC)&alignments_region_mt, 11},
  {"count_region_alignments", (DL_FUNC)&count_region_alignments, 5},
  {"sam_read_n", (DL_FUNC)&sam_read_n, 4},
  {"sam_flag_stats", (DL_FUNC)&sam_flag_stats, 2},
  {"target_lengths", (DL_FUNC)&target_lengths, 1},
  {"bam_cigar_str", (DL_FUNC)&bam_cigar_str, 0},
  {"nuc_table", (DL_FUNC)&nuc_table, 0},
  {"bam_flag", (DL_FUNC)&bam_flag, 1},
  {"extract_aux_tags", (DL_FUNC)&extract_aux_tags, 3},
  {"qpos_to_ref_pos", (DL_FUNC)&qpos_to_ref_pos, 2},
  {"load_bcf", (DL_FUNC)&load_bcf, 2},
  {"bcf_read_n", (DL_FUNC)&bcf_read_n, 3},
  {"arrange_lines", (DL_FUNC)&arrange_lines, 2},
  {NULL, NULL, 0}
};

void R_init_read_bam(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

