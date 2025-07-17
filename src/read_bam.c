#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include "common.h"

// This provides acces to bam / sam and bcf / vcf files.

// The following define what should be returned from calls to
// alignments_region()
// The mess of flags is as usual for historical reasons
#define AR_Q_SEQ 0 // 0x1 // return query seq data
#define AR_Q_DIFF 1 // 0x2 // return positions an identities where query != ref
#define AR_Q_DEPTH 2 // 0x4 // return sequencing depth
#define AR_CIG 3 // 0x8      construct and return a cigar string
#define AR_Q_QUAL 4 // 0x10  return query qualities (ascii encoded; i.e. phred + 0)
#define AR_MT_INFO 5 // 0x20   return information about mate including tlen (not implemented)
#define AR_AUX_MM 6 // 0x40  parse MM info; implies if AR_Q_QUAL, and (AR_Q_SEQ if reference sequence defined).
#define AR_Q_INTRON_DEPTH 7 // 0x80  calculate depths for N operations only; useful for RNA-seq data.

#define MAX_INTRON_L 4096 // if not set by user
// The number of fields in the list returned by alignments_region()
#define AR_R_FIELDS_N 11
static const char* ar_return_fields[AR_R_FIELDS_N] = {"ref", "query", "al", "ops", "seq", "diff", "depth", "cigar", "qual", "mm", "intron.depth"};

// define the bit positions of a set of ret_flag values
// note that these start at 0, so that they can be used
// by: 1 << S_QID  => 0, 1 << S_FLAG => 2, etc..
// This is speicific to sam_read_n() 
#define S_QID 0 // 0x1
#define S_FLAG 1 // 0x2
#define S_RNAME 2  // 0x4
#define S_POS 3  // 0x8
#define S_MAPQ 4 // 0x10
#define S_CIGAR 5 // 0x20
#define S_RNEXT 6 // 0x40
#define S_PNEXT 7 // 0x80
#define S_TLEN 8  // 0x100
#define S_SEQ 9   // 0x200
#define S_QUAL 10  // 0x400
#define S_AUX 11 // 0x800
#define S_CIG_TABLE 12  // 0x1000
#define S_AUX_MM 14  // 0x4000

// to check if a bit is set we can use:
#define bit_set(flag, bit) ( ((1 << (bit)) & (flag)) > 0 )

// and a vector of return types:
//  Note that these are specific to sam_read_n()
// seems reasonable to define the names of these here as well.
// since these depend on each other.
const unsigned int R_ret_types[12] = {STRSXP, INTSXP, INTSXP, INTSXP, INTSXP, // ID, FLAG, RNAME, POS, MAPQ
				      STRSXP, INTSXP, INTSXP, INTSXP, // CIGAR, RNEXT, PNEXT, TLEN
				      STRSXP, STRSXP, STRSXP}; // SEQ, QUAL, AUX


// We want to hold connection data as an external pointer
// to the relevant R data structures:

struct bam_ptrs {
  samFile *sam;
  sam_hdr_t *header;
  hts_idx_t *index;
  // also include an iterator
  hts_itr_t *b_itr;
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

// with functions to read from a bam1_t pointer
struct cigar_string {
  char *buffer;
  size_t string_length; // not including the 0
  size_t buf_size;
};

struct cigar_string cigar_string_init(size_t b_size){
  struct cigar_string cs;
  cs.buffer = malloc(b_size);
  cs.buf_size = b_size;
  cs.string_length = 0;
  return(cs);
}

void cigar_string_free( struct cigar_string *cs ){
  free(cs->buffer);
  cs->buffer = 0;
  cs->buf_size = 0;
  cs->string_length = 0;
}

void cigar_string_grow( struct cigar_string *cs ){
  cs->buf_size = cs->buf_size == 0 ? 256 : cs->buf_size * 2;
  char *nb = malloc(cs->buf_size);
  memcpy( nb, cs->buffer, cs->string_length );
  free(cs->buffer);
  cs->buffer = nb;
}

void cigar_string_read( struct cigar_string *cs, bam1_t *b ){
  if( cs->buf_size == 0 )
    cigar_string_grow(cs);
  cs->string_length = 0;
  cs->buffer[0] = 0; // NULL terminated!
  uint32_t *cig = bam_get_cigar(b);
  int32_t cig_l = b->core.n_cigar;
  for(int i=0; i < cig_l; ++i){
    char op = BAM_CIGAR_STR[ bam_cigar_op(cig[i]) ];
    uint32_t op_l = bam_cigar_oplen( cig[i] );
    while( cs->string_length + 2 + (size_t)log10f((float)op_l) >= cs->buf_size )
      cigar_string_grow( cs );
    cs->string_length += snprintf( cs->buffer + cs->string_length, cs->buf_size - cs->string_length,
				   "%d%c", op_l, op );
  }
}

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
  // It is _not_ safe to: 
  //  free( ptr->sam ); 
  // and there is no hts_file_destroy / sam_file_destroy
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
  warning("Corrupted aux data for read %.*s",
	  b->core.l_qname, bam_get_qname(b));
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
    warning("No more than four modifications per position allowed");
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
// beg:   The first column of the current alignment in the ops table
// end:   The index of the last column of the current alignment + 1
// col_i: The current column being parsed. This should be incremented
//        or decremented depending on the circumstance.
// is_fwd: 0 if the read has been reverse complemented. In that case
//         we have to go backwards.
// qlen:  The length of the query sequence. Needed for reverse mapping.
int query_to_ref(int32_t q_pos, int *r_pos, struct i_matrix *ops,
		 size_t beg, size_t end, size_t *col_i, int is_fwd, int32_t q_len){
  // currently the positions from parse_MM_string are 1 based due
  // to how they are found;
  // however, the positions in the ops table are 0 based.
  // The ops table (since it will be exported to R) should be changed
  // to 1 - based counting at some point in the future. But for now
  // changing all of those is a bit much. Instead here we change
  // the q_pos to a 0 based system.
  //  q_pos--;
  // REMOVE THE ABOVE LINE IF ops table changes to a 1 based system.
  // regression possibility!!
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
  //  int r1 = 5;
  int q1 = 6;
  if(is_fwd){
    while( *col_i < end ){
      int *coords = ops->data + (*col_i) * ops->nrow;
      if( coords[q0] <= q_pos && coords[q1] > q_pos ){
	// if the op type is M, (i.e. 0) we have a match.
	// otherwise we don't have a match.
	if( coords[op] == 0 ){
	  *r_pos = coords[r0] + (q_pos - coords[q0]);
	  return(q_pos);
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
    warning("Obtained a negative query position: %d for q_len %d", q_pos, q_len);
    return(-1);
  }
  while( *col_i >= beg ){
    int *coords = ops->data + (*col_i) * ops->nrow;
    if( coords[q0] <= q_pos && coords[q1] > q_pos ){
      // if the op type is M, (i.e. 0) we have a match.
      // otherwise we don't have a match.
      if( coords[op] == 0 ){
	// do not add 1 here as it's already added to q_pos above.
	*r_pos = coords[r0] + (q_pos - coords[q0]);
	return(q_pos);
      }
      return(q_pos);
    }
    (*col_i)--;
  }
  return(-1);
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
void parse_MM_string(bam1_t *b, uint8_t *s, struct i_matrix *mod_data, int al_i,
		     struct i_matrix *cig_ops, size_t ops_beg, size_t ops_end,
		     const char *ref_seq, size_t ref_seq_l){
  // mod_data should have 6 rows; these are
  // al_i, query_pos, mod_code, mod_n, mod_likelihood, ref_pos
  // This is very wasteful; we could get rid of al_i and return as elements of a list instead
  // but that would just move the problem downstream.
  //
  if(mod_data->nrow != MM_INFO_RN){
    warning("parse_MM_string: mod_data should have %d rows", MM_INFO_RN);
    return;
  }
  if(ops_beg >= ops_end || ops_end > cig_ops->col){
    warning("cigar op rows out of range: %ld -> %ld but only %ld columns of data", ops_beg, ops_end, cig_ops->col);
    return;
  }
  int mod_data_begin = mod_data->col; 
  // return if hard clipped. With hard clipping we have
  // no means of identifying the locations in the sequence.
  uint32_t *cigar = bam_get_cigar(b);
  // Remove any with hard clipping set.
  if(bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP || cig_ops->data[ (ops_end-1) * cig_ops->nrow + 1 ] == BAM_CHARD_CLIP){
    warning("MM parsing not possible for hard clipped sequences");
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
    warning("Modification code data after end of auxiliary string for: %s", bam_get_qname(b));
    return;
  }
  // the type for MM should be Z. We don't know how many positions we are going
  // to parse; 
  // the type for ML should be B, with C, and a 4 bit count 
  // The mm : "MMZ" Followed by:
  // ([ACGTUN][-+]([a-z]+|[0-9]+)[.?]?(,[0-9]+)*;)*
  if(mm[2] != 'Z'){
    warning("Did not find Z after MM");
    return;
  }
  char nuc = mm[3];
  char strand = mm[4];
  // determine how many modifications are related to each entry
  uint32_t mod_code;
  char skip_info;
  int mod_n;
  uint8_t *data = get_ml_code(mm + 5, &mod_code, &mod_n, &skip_info);
  if(data == 0){
    warning("No MM data obtained for query: %s", bam_get_qname(b));
  }
  uint8_t *seq_data = bam_get_seq(b);
  uint8_t *seq_qual = bam_get_qual(b);
  if(seq_data == 0 || seq_qual == 0){
    warning("No sequence or quality data for: %s", bam_get_qname(b));
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
      warning("position exceeded sequence in %s %d --> %d >= %d", bam_get_qname(b), nb, seq_i, seql);
      break;
    }
    // At this point the base at seq_i should be nuc; and should
    // be the one for which we have modification data. 
    // add a column to a table with the following rows:
    // alignment index, query position (seq_i), modification code, mod_n (number of modifications, up to 4), mod_likelihood, ref_pos
    // obtain the refeference position.
    qpos = query_to_ref( seq_i, &r_pos, cig_ops, ops_beg, ops_end, &current_op, is_fwd, seql );
    nuc_info = 0;
    char fwd_base = nuc_encoding[ bam_seqi(seq_data, qpos-1) ];
    if(qpos > 0){
      nuc_info |= (((uint32_t)base) << 24);
      if(seq_qual && qpos <= seql)
	nuc_info |= (((uint32_t)seq_qual[ qpos - 1 ]) << 16);
      if(ref_seq && r_pos <= ref_seq_l){
	nuc_info |= (((uint32_t)ref_seq[ r_pos - 1 ]) << 8);
	nuc_info |= (ref_seq[r_pos-1] == fwd_base);
      }
      // if the base is fwd, set the second bit:
      nuc_info |= (is_fwd << 1);
    }
    push_column( mod_data, (int[]){al_i, (int)qpos, (int)mod_code, mod_n, 0, r_pos, nuc_info} );
    //    push_column( mod_data, (int[]){al_i, (int)seq_i, (int)mod_code, mod_n, 0, r_pos, nuc_info} );
    // find the comma or the semicolon
    while(*data && *data >= '0' && *data <= '9' )
      ++data;
    if(*data == ';'){
      if(data[1] == 0)
	break;
      nuc = data[1];
      strand = data[2];
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
    warning("Expected B in ML auxiliary data");
    return;
  }
  if(ml[3] != 'C'){
    warning("Expected a C in ML auxiliary data but got %c", ml[3]);
    return;
  }
  // An auxiliary with a B type has the following structure:
  // [tag:2 bytes][B: 1 byte][ c|C|s|S|i|I|f 1 byte ][ count: 4 bytes (int) ]
  // Hence we get the number of entries at ml + 4
  uint32_t ml_count = *((int32_t*)(ml + 4));
  if(ml_count != mod_pos_count)
    warning("parse_MM: ml_count (%d) does not equal mod_pos_count (%d) for %s\n(end - ml = %ld)",
	    ml_count, mod_pos_count, bam_get_qname(b), end-ml);
  // And the quality (likelihood of modification data) starts at ml + 8
  uint8_t *qual = ml + 8;
  uint32_t i=0;
  int col = mod_data_begin;
  while(i < ml_count){
    if(col >= mod_data->col){
      warning("exceeded the number of columns in mod_data");
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
  if( col != mod_data->col ){
    warning("col %d is not equal to mod_data->col %ld for %s", col, mod_data->col, bam_get_qname(b));
  }
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
// WARNING: this function counts from 0, not 1; this is different from
//          what I did in aligned_region.
// It is to make it easier to reuse the data within this data.
// Using this, I should be able to subset sequences 
void cigar_to_table(bam1_t *al, int al_i, struct i_matrix *al_coord,
		    int *qcig_length, int *q_beg, int *q_end,
		    size_t *ops_beg, size_t *ops_end){
  // rows are: {"al.i", "op", "type", "r0", "q0", "r1", "q1", "op.l"};
  int column[CIG_OPS_RN] = {0, 0, 0, 0, 0, 0, 0, 0};
  column[0] = al_i;
  uint32_t *cigar = bam_get_cigar(al);
  int32_t cigar_l = al->core.n_cigar;
  // use a 1-based coordinate system instead of a 0-based one
  // for consistency with R and aligned_regions(_
  int r_pos = 1 + (int)al->core.pos;
  int q_pos = 1;
  *qcig_length = 0;
  *q_beg = 1;
  *q_end = 0;
  // if we have hard clipping set q_beg appropriately:
  // we should define MACROs to replace the numbers here.
  if(cigar_l && bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP)
    *q_beg = bam_cigar_oplen(cigar[0]);

  *ops_beg = al_coord->col;
  for(int i=0; i < cigar_l; ++i){
    // if operation is H also increment the qcig_length
    int type = bam_cigar_type(cigar[i]);
    int op = bam_cigar_op(cigar[i]);
    int op_length = bam_cigar_oplen( cigar[i] );
    if( type & 1 || op == BAM_CHARD_CLIP )
      *qcig_length += op_length;
    column[1] = op;
    column[2] = type;
    column[3] = r_pos;
    column[4] = q_pos;
    q_pos += type & 1 ? op_length : 0;
    r_pos += type & 2 ? op_length : 0;
    if(type & 2)
      *q_end = q_pos;
    column[5] = r_pos;
    column[6] = q_pos;
    column[7] = op_length;
    push_column( al_coord, column );
  }
  *ops_end = al_coord->col;
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
//  If ref_seq_r AND (flag_filter & 3 == 3) then return all positions where
//  there is a nucleotide difference
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
  if(region_range[1] < region_range[0] || region_range[0] < 0)
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
  //  if((opt_flag & AR_Q_DIFF) != 0 && (TYPEOF(ref_seq_r) != STRSXP || length(ref_seq_r) != 1))
  if(bit_set(opt_flag, AR_Q_DIFF) && (TYPEOF(ref_seq_r) != STRSXP || length(ref_seq_r) != 1))
    error("Sequence divergence requested but ref_seq is not a character vector of length 1 ");
  const char *ref_seq = (TYPEOF(ref_seq_r) == STRSXP) ? CHAR(STRING_ELT(ref_seq_r, 0)) : 0;
  // which returns an int (which we can change to a size_t)
  size_t ref_seq_l = (size_t)(ref_seq != 0) ? LENGTH(STRING_ELT(ref_seq_r, 0)) : 0;
  //  size_t ref_seq_l = strlen(ref_seq); // this is potentiall rather slow, but necessary
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
  }
  // And similarly for the intron depth:
  int max_intron_length = MAX_INTRON_L; // unless set by user
  if(TYPEOF(max_intron_l_r) == INTSXP && length(max_intron_l_r) == 1){
    max_intron_length = INTEGER(max_intron_l_r)[0];
  }else{
    warning("max_intron_length set to default value (%d)", max_intron_length);
  }
  int *intron_depth = 0;
  if(bit_set(opt_flag, AR_Q_INTRON_DEPTH)){
    SET_VECTOR_ELT( ret_data, 10, allocVector(INTSXP, region_length));
    intron_depth = INTEGER(VECTOR_ELT(ret_data, 10));
    memset( intron_depth, 0, sizeof(int) * region_length );
  }
  // We use sam_itr_queryi, so that we can specify the region
  // numerically thus also allowing the same function to return the sequence depth for the region.
  hts_itr_t *b_itr = sam_itr_queryi(bam->index, target_id, region_range[0], region_range[1]);
  // hts_itr_destroy( ) not called on this anywhere: suggests we have a memory leak.
  // previously I used the string version:
  //  hts_itr_t *b_itr = sam_itr_querys(bam->index, bam->header, region);
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
  
  //  const char* rownames[8] = {"al.i", "op", "type", "r0", "q0", "r1", "q1", "op.l"};
  // struct i_matrix al_coord = init_i_matrix( 8, init_size );
  // rownames are in cig_ops_rownames
  struct i_matrix al_coord = init_i_matrix( CIG_OPS_RN, init_size );
  SEXP al_coord_rownames_r = PROTECT(mk_rownames( cig_ops_rownames, CIG_OPS_RN ));

  // For sequence variants in the reads:
  //  const char* var_coord_row_names[] = {"seq.id", "r.pos", "q.pos", "nuc"};
  struct i_matrix var_coord = init_i_matrix( 4, init_size );
  SEXP var_coord_row_names_r = PROTECT(mk_rownames( (const char*[]){"al.i", "r.pos", "q.pos", "nuc"}, 4 ));

  // For base modification data:
  struct i_matrix mm_info = init_i_matrix( MM_INFO_RN, init_size );
  SEXP mm_info_rownames_r = PROTECT(mk_rownames( mm_info_rownames, MM_INFO_RN ));
  
  //  int column[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  // dimNames_r used to set rownames for the main return table
  //  SEXP dimNames_r = PROTECT(mk_rownames( rownames, 8 ));
  // We also want to return more basic information about the alignments eg.
  // flag, pos, mapq, qlen and tlen
  struct i_matrix al_vars = init_i_matrix( 8, init_size );
  int av_column[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  SEXP al_var_dimnames = PROTECT( mk_rownames((const char*[]){"flag", "r.beg", "r.end", "q.beg", "q.end", "mqual", "qlen", "qclen"}, 8));
  char *qseq = 0; // use as a temporary variable to hold the address of the query sequence
  // Like qseq, but for the quality values. These can be obtained using
  // bam_get_qual() which returns a pointer to the qual data
  //                of length bam->core.l_qseq
  // Currently there is no option to return the quality scores; but I can add qualities
  // at mismatched position without modifying any data structures. I
  unsigned char *qqual=0; //
  while((r=sam_itr_next(bam->sam, b_itr, al)) >= 0){
    qseq = qqual = 0;
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
	str_array_push_cp_n(&query_qual, qqual, al->core.l_qseq);
      else
	str_array_push_cp(&query_qual, "*");
    }
    al_count++;
    //    column[0] = al_count;
    // then go through the cigar string and for each op create a new column
    // of offsets for the al_count
    uint32_t *cigar = bam_get_cigar(al);
    int32_t cigar_l = al->core.n_cigar;
    int r_pos = 1 + (int)al->core.pos;
    int q_pos = 1;
    // We don't actually need the av_column; we could simply use an anonymous
    // array.
    av_column[0] = (int)al->core.flag;
    av_column[1] = r_pos;
    // 2, 3, 4 are r.end, q.beg, q.end respectively
    av_column[5] = (int)al->core.qual;
    av_column[6] = (int)al->core.l_qseq;
    // It's not actually clear what isize is? Seems to be set to 0
    // a lot of the time.. set to qcig_length instead
    // av_column[7] = (int)al->core.isize;
    int q_beg=-1;
    int q_end=0;
    // calculate the true query length from the cigar string
    int qcig_length = 0;
    // In order to parse MM data we need to pass the columns of the al_coord table
    // that refer to the current alignment;
    size_t ops_begin_col = al_coord.col;
    for(int i=0; i < cigar_l; ++i){
      // if operation is H also increment the qcig_length
      if( bam_cigar_type(cigar[i]) & 1 || bam_cigar_op( cigar[i] ) == BAM_CHARD_CLIP )
	qcig_length += bam_cigar_oplen( cigar[i] );
      int cig_op = bam_cigar_op(cigar[i]);
      int cig_type = bam_cigar_type(cigar[i]);
      int cig_oplen = bam_cigar_oplen(cigar[i]);
      int ref_0 = r_pos;
      int ref_1 = r_pos + (cig_type & CIG_OP_TP_RC ? cig_oplen : 0);
      int q_0 = q_pos;
      int q_1 = q_pos +  (cig_type & CIG_OP_TP_QC ? cig_oplen : 0);
      if( q_beg == -1 && (cig_type & CIG_OP_TP_RC) )
	q_beg = q_pos;
      // If the caller has requested the variance information and the
      // cigar op consumes both the query and the reference
      // then check for differences in the sequence
      if( bit_set(opt_flag, AR_Q_DIFF) || bit_set(opt_flag, AR_Q_DEPTH) && cig_type == 3 ){
	// I added +1 to r_pos and q_pos above; this means that I
	// need to subtract 1 from all of j here.
	for(int j=0; j < cig_oplen; ++j){
	  // Note that qseq may not be stored.
	  if(qseq && (bit_set(opt_flag, AR_Q_DIFF)) && (j + r_pos) < ref_seq_l ){
	    if(qseq[ q_pos + j - 1] != ref_seq[ r_pos + j -1 ]){
	      int nuc_info = ((int)ref_seq[r_pos + j -1 ] << 8) | (int)qseq[q_pos+j -1];
	      if(qqual)
		nuc_info |= (qqual[q_pos+j-1] << 16);
	      // This is a bit of a nameful, but see definition of data struct above
	      push_column(&var_coord, (int[]){al_count, r_pos+j, q_pos+j,
						nuc_info});
	    }
	  }
	  if(bit_set(opt_flag, AR_Q_DEPTH)){
	    // r_pos is 1 based
	    int o = (r_pos+j-1) - region_range[0];
	    if(o >= 0 && o < region_length)
	      seq_depth[o]++;
	  }
	}
      }
      if( (bit_set(opt_flag, AR_Q_INTRON_DEPTH)) && cig_op == BAM_CREF_SKIP
	  && cig_oplen <= max_intron_length){
	for(int j=0; j < cig_oplen; ++j){
	  // r_pos is 1 based
	  int o = (r_pos+j-1) - region_range[0];
	  if(o >= 0 && o < region_length)
	    intron_depth[o]++;
	}
      }
      // if an insertion then we need to increment the depth at a single
      // reference position
      // WRONG; this is not needed because the next operation will be a match
      // op starting at the same position. This means that we will overcount the depth.
      q_pos += bam_cigar_type(cigar[i]) & CIG_OP_TP_QC ? bam_cigar_oplen( cigar[i] ) : 0;
      r_pos += bam_cigar_type(cigar[i]) & CIG_OP_TP_RC ? bam_cigar_oplen( cigar[i] ) : 0;
      if(bam_cigar_type(cigar[i]) & 2)
	q_end = q_pos;
      push_column( &al_coord, (int[]){al_count, cig_op, cig_type, ref_0, q_0, ref_1, q_1, cig_oplen} );
    }
    size_t ops_end_col = al_coord.col;
    av_column[2] = r_pos;
    av_column[3] = q_beg;
    av_column[4] = q_end;
    av_column[7] = qcig_length;
    push_column(&al_vars, av_column);
    // We should modify the function so that we do not malloc and free for every query
    // sequence.
    if(bit_set(opt_flag, AR_AUX_MM)){
      parse_MM_string(al, 0, &mm_info, al_count, &al_coord, ops_begin_col, ops_end_col, ref_seq, ref_seq_l);
    }
    // qseq may be 0, but that should be safe to free.
    // We only need to free qseq if we did not store it in the str_array_object
    // In fact we should never need to do this.
    if(!bit_set(opt_flag, AR_Q_SEQ) )
      free(qseq);
    qseq = 0;
  }
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
  /* if(TYPEOF(min_ql_r) != INTSXP || length(min_ql_r) != 1) */
  /*   error("min_ql_r should be an integer vector of length 1"); */
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
  // To start with to make sure that the function works, simply return the identifiers.
  // afterwards we can consider more complex options. But for this we mostly will want
  // the id, flag, cigar, seq and quality.. But we can check that later.
  //  SEXP ret_data = PROTECT( allocVector(VECSXP, 5) );
  // n, the number of reads obtained, and then the 12 fields of the sam forma;
  SEXP ret_data_names = PROTECT( mk_strsxp( (const char*[]){"id", "flag", "ref", "pos", "mapq",
							      "cigar", "ref.m", "pos.m", "tlen",
							      "seq", "qual", "aux", "ops", "q.inf", "mm", "n"}, 16));
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
      cigar_to_table(b, 1+(*count), &al_coord, &qcig_length, &q_beg, &q_end, &ops_begin, &ops_end);
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
  {"alignments_region", (DL_FUNC)&alignments_region, 9},
  {"count_region_alignments", (DL_FUNC)&count_region_alignments, 5},
  {"sam_read_n", (DL_FUNC)&sam_read_n, 4},
  {"sam_flag_stats", (DL_FUNC)&sam_flag_stats, 2},
  {"target_lengths", (DL_FUNC)&target_lengths, 1},
  {"bam_cigar_str", (DL_FUNC)&bam_cigar_str, 0},
  {"nuc_table", (DL_FUNC)&nuc_table, 0},
  {"bam_flag", (DL_FUNC)&bam_flag, 1},
  {"extract_aux_tags", (DL_FUNC)&extract_aux_tags, 3},
  {"load_bcf", (DL_FUNC)&load_bcf, 2},
  {"bcf_read_n", (DL_FUNC)&bcf_read_n, 3},
  {"arrange_lines", (DL_FUNC)&arrange_lines, 2},
  {NULL, NULL, 0}
};

void R_init_read_bam(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

