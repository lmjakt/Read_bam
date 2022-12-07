#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <htslib/sam.h>

// not sure why this isn't in a header somewhere
// but I can't find it.
// this is a copy of seq_nt16_str defined as a non null-
// terminated vector array
const char *nuc_encoding = "=ACMGRSVTWYHKDBN";


// We want to hold connection data as an external pointer
// to the relevant R data structures:

struct bam_ptrs {
  samFile *sam;
  sam_hdr_t *header;
  hts_idx_t *index;
};

// that will be held as an external pointer
// Which should be cleaned up

static void finalise_bam_ptrs(SEXP ptr_r){
  //  Rprintf("finalising bam_ptrs\n");
  if(!R_ExternalPtrAddr(ptr_r))  // if already NULL, do nothing
    return;
  struct bam_ptrs *ptr = (struct bam_ptrs*)R_ExternalPtrAddr(ptr_r);
  // destroy the various indices
  sam_hdr_destroy( ptr->header );
  hts_idx_destroy( ptr->index );
  // It is _not_ safe to: 
  //  free( ptr->sam ); 
  // and there is no hts_file_destroy / sam_file_destroy
  free( ptr ); // ??
  R_ClearExternalPtr(ptr_r);
}

// These are structures for dynamically growing a return data set;
// They would be better off in a library somewhere as I keep repeating
// this code.

// row major? as in R
struct i_matrix {
  int *data;
  size_t nrow;
  size_t ncol;
  size_t row;
  size_t col;
};

struct i_matrix init_i_matrix(size_t nrow, size_t ncol){
  struct i_matrix m;
  m.nrow=nrow;
  m.ncol=ncol;
  m.row=0;
  m.col=0;
  m.data = malloc( sizeof(int) * m.nrow * m.ncol );
  return(m);
}

void double_columns( struct i_matrix *m ){
  int *old = m->data;
  size_t sz = sizeof(int) * m->nrow * m->ncol;
  m->data = malloc( 2 * sz );
  memcpy( m->data, old, sz );
  m->ncol = m->ncol * 2;
  free(old);
}

// push numbers in c onto the matrix
void push_column( struct i_matrix *m, int *c ){
  if(m->col >= m->ncol)
    double_columns( m );
  memcpy( m->data + m->col * m->nrow, c, sizeof(int) * m->nrow );
  m->col++;
}


struct str_array {
  size_t capacity;
  size_t length;
  char **strings;
};

struct str_array init_str_array(size_t init_size){
  struct str_array str;
  str.capacity = init_size;
  str.length = 0;
  str.strings = malloc(sizeof(char*) * init_size);
  return(str);
}

void double_str_array(struct str_array *str){
  char **old = str->strings;
  str->capacity = str->capacity * 2;
  str->strings = malloc(sizeof(char*) * str->capacity );
  memcpy( str->strings, old, sizeof(char*) * str->length );
  free(old);
}

void str_array_push_cp(struct str_array *str, const char *word){
  if(str->length >= str->capacity)
    double_str_array(str);
  size_t l = strlen(word);
  str->strings[str->length] = malloc(sizeof(char) * (l + 1));
  memcpy( str->strings[str->length], word, sizeof(char) * (l + 1));
  str->length++;
}

void str_array_push(struct str_array *str, char *word){
  if(str->length >= str->capacity)
    double_str_array(str);
  str->strings[str->length] = word;
  str->length++;
}


void str_array_free(struct str_array *str){
  for(size_t i=0; i < str->length; ++i)
    free(str->strings[i]);
  free(str->strings);
}

// Convenience function that makes a char sequence out of
// the sequence in bam1_t.data

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
    error("Unable to open index");
  }
  sam_hdr_t *header = sam_hdr_read(sam);
  // Do I need to allocate this on the stack?
  struct bam_ptrs *ptr = malloc( sizeof(struct bam_ptrs) );
  ptr->sam = sam;
  ptr->header = header;
  ptr->index = index;
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

struct bam_ptrs* extract_bam_ptr(SEXP bam_ptr_r){
  if(TYPEOF(bam_ptr_r) != EXTPTRSXP)
    error("bam_ptr_r should be an external pointer");
  SEXP bam_tag = PROTECT(R_ExternalPtrTag(bam_ptr_r));
  //  SEXP bam_prot = PROTECT(R_ExternalPtrProtected(bam_ptr_r));
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
//  If ref_seq_r AND (flag_filter & 3 == 3) then return all positions where
//  there is a nucleotide difference
//  min_mq_r : A minimum mapping quality
//  min_ql_r : A minimum query length, taken from bam1_core_t.qlen
//             which does not necessarily equate to the read length
SEXP alignments_region(SEXP region_r, SEXP region_range_r,
		       SEXP bam_ptr_r, SEXP flag_filter_r, 
		       SEXP opt_flag_r, SEXP ref_seq_r,
		       SEXP min_mq_r, SEXP min_ql_r){
  if(TYPEOF(region_r) != STRSXP || length(region_r) < 1)
    error("region_r should be a character vector of positive length");
  struct bam_ptrs *bam = extract_bam_ptr(bam_ptr_r);
  if(!bam){
    // we could be clever here and try to recreate from the protect information
    // but we should think carefully about this
    error("External pointer is NULL");
  }
  if(TYPEOF(region_range_r) != INTSXP || length(region_range_r) != 2)
    error("region_range should be an integer vector of length 2");
  int *region_range = INTEGER(region_range_r);
  if(region_range[1] < region_range[0] || region_range[0] < 0)
    error("Unsorted region_range");
  const char *region = CHAR(STRING_ELT(region_r, 0));
  int target_id = sam_hdr_name2tid(bam->header, region);
  if(target_id < 0)
    error("Unable to find target id for specified region: %s", region);
  
  if(TYPEOF(flag_filter_r) != INTSXP || length(flag_filter_r) != 2)
    error("flag_filter should be an integer vector of length 2");
  if(TYPEOF(opt_flag_r) != INTSXP || length(opt_flag_r) != 1)
    error("opt_flag should be a single integer value");
  int opt_flag = asInteger(opt_flag_r);
  if(opt_flag < 0)
    opt_flag = 0;
  if((opt_flag & 2) == 2 && (TYPEOF(ref_seq_r) != STRSXP || length(ref_seq_r) != 1))
    error("Sequence divergence requested but ref_seq is not a character vector of length 1 ");
  const char *ref_seq = ((opt_flag & 2) == 2) ? CHAR(STRING_ELT(ref_seq_r, 0)) : "";
  size_t ref_seq_l = strlen(ref_seq); // this is potentiall rather slow, but necessary
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
  SEXP ret_data = PROTECT(allocVector(VECSXP, 7));
  SEXP ret_names_r = PROTECT( mk_strsxp( (const char*[]){"ref", "query", "al", "ops", "seq", "diff", "depth"}, length(ret_data) ));
  setAttrib( ret_data, R_NamesSymbol, ret_names_r );

  // If depth requested we can add this to the data structure immediately..
  int *seq_depth = 0;
  int region_length = 1 + region_range[1] - region_range[0];
  if(opt_flag & 4){
    SET_VECTOR_ELT( ret_data, 6, allocVector(INTSXP, region_length) );
    seq_depth = INTEGER(VECTOR_ELT(ret_data, 6));
    memset( seq_depth, 0, sizeof(int) * region_length );
  }
  
  // and then we do the whole business of getting the reads..
  // We want to change this to use sam_itr_queryi, so that we specify the region
  // numerically thus also allowing the same function to return the sequence depth for the region.
  hts_itr_t *b_itr = sam_itr_queryi(bam->index, target_id, region_range[0], region_range[1]);
  // previously I used the string version:
  //  hts_itr_t *b_itr = sam_itr_querys(bam->index, bam->header, region);
  bam1_t *al = bam_init1();
  int r=0;
  int al_count = 0;
  // Data structures that will hold the return data:
  size_t init_size = 10;
  struct str_array query_ids = init_str_array(init_size);
  // query_seq is filled ony if opt_flag & 1 == 1
  struct str_array query_seq = init_str_array(init_size);
  const char* rownames[7] = {"al.i", "op", "type", "r0", "q0", "r1", "q1"};
  struct i_matrix al_coord = init_i_matrix( 7, init_size );

  // For sequence variants in the reads:
  //  const char* var_coord_row_names[] = {"seq.id", "r.pos", "q.pos", "nuc"};
  struct i_matrix var_coord = init_i_matrix( 4, init_size );
  SEXP var_coord_row_names_r = PROTECT(mk_rownames( (const char*[]){"seq.id", "r.pos", "q.pos", "nuc"}, 4 ));

  int column[7] = {0, 0, 0, 0, 0, 0, 0};
  // dimNames_r used to set rownames for the main return table
  SEXP dimNames_r = PROTECT(mk_rownames( rownames, 7 ));
  // We also want to return more basic information about the alignments eg.
  // flag, pos, mapq, qlen and tlen
  struct i_matrix al_vars = init_i_matrix( 8, init_size );
  int av_column[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  SEXP al_var_dimnames = PROTECT( mk_rownames((const char*[]){"flag", "r.beg", "r.end", "q.beg", "q.end", "mqual", "qlen", "tlen"}, 8));
  char *qseq = 0; // use as a temporary variable to hold the address of the query sequence
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
    if(al->core.l_qseq < min_ql)
      continue;
    str_array_push_cp(&query_ids, bam_get_qname(al));
    if(opt_flag & 3){
      qseq = bam_seq(al);
      if(opt_flag & 1){
	if(qseq)
	  str_array_push(&query_seq, qseq);
	else
	  str_array_push_cp(&query_seq, "*");
      }
    }
    al_count++;
    column[0] = al_count;
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
    av_column[7] = (int)al->core.isize;
    int q_beg=-1;
    int q_end=0;
    for(int i=0; i < cigar_l; ++i){
      column[1] = bam_cigar_op( cigar[i] ) + 1;
      column[2] = bam_cigar_type(cigar[i]);
      column[3] = r_pos;
      column[4] = q_pos;
      if( q_beg == -1 && (bam_cigar_type(cigar[i]) & 2) )
	q_beg = q_pos;
      // If the caller has requested the variance information and the
      // cigar op consumes both the query and the reference
      // then check for differences in the sequence
      if( (opt_flag & 6) &&  bam_cigar_type(cigar[i]) == 3){
	// I added +1 to r_pos and q_pos above; this means that I
	// need to subtract 1 from all of j here.
	for(int j=0; j < bam_cigar_oplen(cigar[i]); ++j){
	  // Note that qseq may not be stored.
	  if(qseq && (opt_flag & 2) && (j + r_pos) < ref_seq_l ){
	    if(qseq[ q_pos + j - 1] != ref_seq[ r_pos + j -1 ]){
	      // This is a bit of a nameful, but see definition of data struct above
	      push_column(&var_coord, (int[]){al_count, r_pos+j, q_pos+j,
		    (((int)ref_seq[r_pos + j -1 ] << 8) | (int)qseq[q_pos+j -1])} );
	    }
	  }
	  if(opt_flag & 4){
	    int o = (r_pos+j) - region_range[0];
	    if(o >= 0 && o < region_length)
	      seq_depth[o]++;
	  }
	}
      }
      q_pos += bam_cigar_type(cigar[i]) & 1 ? bam_cigar_oplen( cigar[i] ) : 0;
      r_pos += bam_cigar_type(cigar[i]) & 2 ? bam_cigar_oplen( cigar[i] ) : 0;
      if(bam_cigar_type(cigar[i]) & 2)
	q_end = q_pos;
      column[5] = r_pos;
      column[6] = q_pos;
      push_column( &al_coord, column );
    }
    av_column[2] = r_pos;
    av_column[3] = q_beg;
    av_column[4] = q_end;
    push_column(&al_vars, av_column);
    // qseq may be 0, but that should be safe to free.
    // We only need to free qseq if we did not store it in the str_array_object
    if( (opt_flag & 1) == 0 )
      free(qseq);
    qseq = 0;
  }
  // query ids
  SET_VECTOR_ELT(ret_data, 0, region_r);
  SET_VECTOR_ELT(ret_data, 1, allocVector(STRSXP, query_ids.length));
  // Handle query ids and query sequences at the same time as they have
  // the same length.
  if(opt_flag & 1)
    SET_VECTOR_ELT(ret_data, 4, allocVector(STRSXP, query_ids.length));
  SEXP q_ids_r = VECTOR_ELT( ret_data, 1 );
  for(size_t i=0; i < query_ids.length; ++i){
    SET_STRING_ELT( q_ids_r, i, mkChar( query_ids.strings[i] ));
    if(opt_flag & 1)
      SET_STRING_ELT( VECTOR_ELT(ret_data, 4), i, mkChar( query_seq.strings[i] ));
  }
  str_array_free(&query_ids);
  str_array_free(&query_seq);
  // alignment variables
  SET_VECTOR_ELT(ret_data, 2, allocMatrix(INTSXP, al_vars.nrow, al_vars.col));
  setAttrib( VECTOR_ELT(ret_data, 2), R_DimNamesSymbol, al_var_dimnames );
  memcpy( INTEGER(VECTOR_ELT(ret_data, 2)), al_vars.data, sizeof(int) * al_vars.nrow * al_vars.col );
  free(al_vars.data);
  // alignment coordinates
  SET_VECTOR_ELT(ret_data, 3, allocMatrix(INTSXP, al_coord.nrow, al_coord.col));
  //  setAttrib( VECTOR_ELT(ret_data, 1), R_RowNamesSymbol, rowNames_r);
  setAttrib( VECTOR_ELT(ret_data, 3), R_DimNamesSymbol, dimNames_r);
  memcpy( INTEGER(VECTOR_ELT(ret_data, 3)), al_coord.data, sizeof(int) * al_coord.nrow * al_coord.col);
  free(al_coord.data);
  // If we have variant data then we need to do something with it before freeing up the used memory.
  if(opt_flag & 2){
    SET_VECTOR_ELT(ret_data, 5, allocMatrix(INTSXP, var_coord.nrow, var_coord.col));
    setAttrib( VECTOR_ELT(ret_data, 5), R_DimNamesSymbol, var_coord_row_names_r );
    memcpy( INTEGER(VECTOR_ELT(ret_data, 5)), var_coord.data, sizeof(int) * var_coord.nrow * var_coord.col );
  }
  free(var_coord.data);
  UNPROTECT(5);
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


static const R_CallMethodDef callMethods[] = {
  {"load_bam", (DL_FUNC)&load_bam, 2},
  {"alignments_region", (DL_FUNC)&alignments_region, 8},
  {"target_lengths", (DL_FUNC)&target_lengths, 1},
  {"bam_cigar_str", (DL_FUNC)&bam_cigar_str, 0},
  {"nuc_table", (DL_FUNC)&nuc_table, 0},
  {"bam_flag", (DL_FUNC)&bam_flag, 1},
  {NULL, NULL, 0}
};

void R_init_read_bam(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

