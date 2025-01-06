#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <htslib/sam.h>
#include "split_string.h"

// define the bit positions of a set of ret_flag values
// note that these start at 0, so that they can be used
// by: 1 << S_QID  => 0, 1 << S_FLAG => 2, etc..
#define S_QID 0
#define S_FLAG 1
#define S_RNAME 2
#define S_POS 3
#define S_MAPQ 4
#define S_CIGAR 5
#define S_RNEXT 6
#define S_PNEXT 7
#define S_TLEN 8
#define S_SEQ 9
#define S_QUAL 10
#define S_AUX 11
#define S_CIG_TABLE 12

// The number of rows in an cigar operations table
#define CIG_OPS_RN 8
#define OPS_INIT_SIZE 256
// a global variable; I feel dirty;
const char* cig_ops_rownames[CIG_OPS_RN] = {"al.i", "op", "type", "r0", "q0", "r1", "q1", "op.l"};

// the following is not that useful.. 
// #define BIT_F( F ) ( 1 << F )

// and a vector of return types:
const unsigned int R_ret_types[12] = {STRSXP, INTSXP, INTSXP, INTSXP, INTSXP, // ID, FLAG, RNAME, POS, MAPQ
				      STRSXP, INTSXP, INTSXP, INTSXP, // CIGAR, RNEXT, PNEXT, TLEN
				      STRSXP, STRSXP, STRSXP}; // SEQ, QUAL, AUX

// not sure why this isn't in a header somewhere
// but I can't find it.
// this is a copy of seq_nt16_str defined as a non null-
// terminated vector array
const char *nuc_encoding = "=ACMGRSVTWYHKDBN";
// if we want to reverse complement
const char *nuc_encoding_rc = "=TGKCYSBAWRDMHVN";

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
  if(m->ncol == 0){
    m->ncol = 2;
    m->data = malloc(sizeof(int) * m->nrow * m->ncol);
    return;
  }
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

// Avoid malloc if possible:
// but note that whatever we do, we end up copying data too much
// Note that we could probably use, a single call of `nibble2base`
// to get the sequence.
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

// Parse the cigar data and extend:
// A cigar operations table containing the operations, their lengths and
//     the query and reference positions (i.e. how they align)
// Sets q_length the query length calculated from the cigar string; this can differ
//     from that given in the bam1_t data struct due to hard clipping.
//     Similarly, the beginning and end positions of the original query can differ
//     So the function also sets q_beg and q_end. These refer to the query
//     prior to any hard clipping.
// 
// In the future I might add functionality that also identifies sequence
// differences and parses MM and ML data; however, it is difficult to do this
// without making spaghetti code;
//
// WARNING: this function counts from 0, not 1; this is different from
//          what I did in aligned_region.
// It is to make it easier to reuse the data within this data.
// Using this, I should be able to subset sequences 
void cigar_to_table(bam1_t *al, int al_i, struct i_matrix *al_coord, int *qcig_length, int *q_beg, int *q_end){
  // rows are: {"al.i", "op", "type", "r0", "q0", "r1", "q1", "op.l"};
  int column[CIG_OPS_RN] = {0, 0, 0, 0, 0, 0, 0, 0};
  column[0] = al_i;
  uint32_t *cigar = bam_get_cigar(al);
  int32_t cigar_l = al->core.n_cigar;
  int r_pos = (int)al->core.pos;
  int q_pos = 0;
  *qcig_length = 0;
  *q_beg = 1;
  *q_end = 0;
  // if we have hard clipping set q_beg appropriately:
  if(cigar_l && bam_cigar_op(cigar[0]) == 5)
    *q_beg = bam_cigar_oplen(cigar[0]);
  
  for(int i=0; i < cigar_l; ++i){
    // if operation is H also increment the qcig_length
    int type = bam_cigar_type(cigar[i]);
    int op = bam_cigar_op(cigar[i]);
    int op_length = bam_cigar_oplen( cigar[i] );
    if( type & 1 || op == 5 )
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
  SEXP ret_names_r = PROTECT( mk_strsxp( (const char*[]){"ref", "query", "al", "ops", "seq", "diff", "depth", "cigar"}, 8 ));
  SEXP ret_data = PROTECT(allocVector(VECSXP, length(ret_names_r)));
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
  size_t init_size = 10;
  struct str_array query_ids = init_str_array(init_size);
  // query_seq is filled ony if opt_flag & 1 == 1
  struct str_array query_seq = init_str_array(init_size);
  // We unfortunaly have ended up with a struct cigar_string,
  // and cigar, cigars and cigars_string variables; That is BAD, and ought to be cleaned up
  struct str_array cigars = init_str_array(init_size);
  struct cigar_string cigars_string = cigar_string_init(256); // I feel dirty..
  const char* rownames[8] = {"al.i", "op", "type", "r0", "q0", "r1", "q1", "op.l"};
  struct i_matrix al_coord = init_i_matrix( 8, init_size );

  // For sequence variants in the reads:
  //  const char* var_coord_row_names[] = {"seq.id", "r.pos", "q.pos", "nuc"};
  struct i_matrix var_coord = init_i_matrix( 4, init_size );
  SEXP var_coord_row_names_r = PROTECT(mk_rownames( (const char*[]){"seq.id", "r.pos", "q.pos", "nuc"}, 4 ));

  int column[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  // dimNames_r used to set rownames for the main return table
  SEXP dimNames_r = PROTECT(mk_rownames( rownames, 8 ));
  // We also want to return more basic information about the alignments eg.
  // flag, pos, mapq, qlen and tlen
  struct i_matrix al_vars = init_i_matrix( 8, init_size );
  int av_column[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  SEXP al_var_dimnames = PROTECT( mk_rownames((const char*[]){"flag", "r.beg", "r.end", "q.beg", "q.end", "mqual", "qlen", "qclen"}, 8));
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
    // read the cigar and push it to cigars
    if(opt_flag & 8){
      cigar_string_read( &cigars_string, al );
      str_array_push_cp(&cigars, cigars_string.buffer );
    }
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
    // It's not actually clear what isize is? Seems to be set to 0
    // a lot of the time.. set to qcig_length instead
    // av_column[7] = (int)al->core.isize;
    int q_beg=-1;
    int q_end=0;
    // calculate the true query length from the cigar string
    int qcig_length = 0;
    for(int i=0; i < cigar_l; ++i){
      // if operation is H also increment the qcig_length
      if( bam_cigar_type(cigar[i]) & 1 || bam_cigar_op( cigar[i] ) == 5 )
	qcig_length += bam_cigar_oplen( cigar[i] );
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
      // if an insertion then we need to increment the depth at a single
      // reference position
      // WRONG; this is not needed because the next operation will be a match
      // op starting at the same position. This means that we will overcount the depth.
      q_pos += bam_cigar_type(cigar[i]) & 1 ? bam_cigar_oplen( cigar[i] ) : 0;
      r_pos += bam_cigar_type(cigar[i]) & 2 ? bam_cigar_oplen( cigar[i] ) : 0;
      if(bam_cigar_type(cigar[i]) & 2)
	q_end = q_pos;
      column[5] = r_pos;
      column[6] = q_pos;
      column[7] = bam_cigar_oplen( cigar[i] );
      push_column( &al_coord, column );
    }
    av_column[2] = r_pos;
    av_column[3] = q_beg;
    av_column[4] = q_end;
    av_column[7] = qcig_length;
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
  if(opt_flag & 8)
    SET_VECTOR_ELT(ret_data, 7, allocVector(STRSXP, cigars.length));
  // Handle query ids and query sequences at the same time as they have
  // the same length.
  if(opt_flag & 1)
    SET_VECTOR_ELT(ret_data, 4, allocVector(STRSXP, query_ids.length));
  SEXP q_ids_r = VECTOR_ELT( ret_data, 1 );
  SEXP cigars_r = VECTOR_ELT( ret_data, 7 );
  for(size_t i=0; i < query_ids.length; ++i){
    SET_STRING_ELT( q_ids_r, i, mkChar( query_ids.strings[i] ));
    if(opt_flag & 8)
      SET_STRING_ELT( cigars_r, i, mkChar( cigars.strings[i] ));
    if(opt_flag & 1)
      SET_STRING_ELT( VECTOR_ELT(ret_data, 4), i, mkChar( query_seq.strings[i] ));
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
							      "seq", "qual", "aux", "ops", "q.inf", "n"}, 15));
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
  // The last element may be a matrix, if specified. Otherwise left as NULL.
  // But we have to copy the matrix as it contains all the cigar operations
  // and these are of unknown length;
  // see top of the file for what these are.. 
  // The following four of these are character vectors: id, seq, qual, aux
  /* for(int i=1; i < length(ret_data); ++i) */
  /*   SET_VECTOR_ELT(ret_data, i, allocVector(STRSXP, n)); */
  //  int r;  // the return value from sam_read1
  bam1_t *b = bam_init1();
  size_t seq_buffer_size = 500;
  char *seq_buffer = malloc(seq_buffer_size);
  kstring_t aux_str = KS_INITIALIZE;
  // If S_CIG_TABLE is set, then initialise to i_matrix arrays;
  // one will hold q_info (the actual query lengths, beg and end positions before clipping)
  // and the other will hold the cigar operations in a table of unknown length;
  int ret_cig_ops = ret_flag & (1 << S_CIG_TABLE);
  struct i_matrix al_coord = init_i_matrix(CIG_OPS_RN, (ret_cig_ops ? OPS_INIT_SIZE : 0));
  // they query_info will hold the q_length, q_cig_length, q_beg, q_end. (possibly also
  struct i_matrix query_info = init_i_matrix(4, (ret_cig_ops ? n : 0));
  // We should then free the al_coord.data query_info.data after copying to an R data
  // structure.
  // the initial size should probably be exposed to the user
  struct cigar_string cigar = cigar_string_init(256);
  // size_t aux_str_size = 0;  return value not used.
  //  for(int i=0; i < n && sam_read1(bam->sam, bam->header, b) >= 0; ++i){
  int qcig_length, q_beg, q_end;
  while( (*count) < n && read_next_entry(bam, b) >= 0){
    //  for(int i=0; i < n && read_next_entry(bam, b) >= 0; ++i){
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
      cigar_to_table(b, 1+(*count), &al_coord, &qcig_length, &q_beg, &q_end);
      push_column(&query_info, (int[]){ b->core.l_qseq, qcig_length, q_beg, q_end });
    }
    (*count)++;
  }
  // if S_CIG_TABLE, make the appropriate data structures:
  if(ret_flag & (1 << S_CIG_TABLE)){
    SET_VECTOR_ELT( ret_data, S_CIG_TABLE, allocMatrix(INTSXP, al_coord.nrow, al_coord.col) );
    setAttrib(VECTOR_ELT(ret_data, S_CIG_TABLE), R_DimNamesSymbol, mk_rownames(cig_ops_rownames, CIG_OPS_RN));
    SET_VECTOR_ELT( ret_data, S_CIG_TABLE+1, allocMatrix(INTSXP, query_info.nrow, query_info.col) );
    setAttrib(VECTOR_ELT(ret_data, S_CIG_TABLE+1), R_DimNamesSymbol,
	      mk_rownames((const char*[]){"qlen", "q.cigl", "q.beg", "q.end"}, 4));
    memcpy( INTEGER(VECTOR_ELT(ret_data, S_CIG_TABLE)), al_coord.data, sizeof(int) * al_coord.nrow * al_coord.col );
    memcpy( INTEGER(VECTOR_ELT(ret_data, S_CIG_TABLE+1)), query_info.data, sizeof(int) * query_info.nrow * query_info.col );
  }
  free(seq_buffer);
  free(aux_str.s);
  free(al_coord.data);
  free(query_info.data);
  cigar_string_free( &cigar );
  bam_destroy1(b);
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
  if(TYPEOF(aux_types_r) != STRSXP || length(aux_types_r) < 1)
    error("aux_types_r must be a string with length of 1 or longer");
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
  // which give values of, 3, 19, 9, 6, 26, and 8 respectively.
  // That means we can check if a type is defined by AND operations on a 32 bit unsigned
  // integer
  // note that toupper is equivalent to & 0xDF
  uint8_t up_mask = 0xDF;
  // This is with bits 17, 3, 19, 9, 6, 26 and 8 set.
  uint32_t type_check = 0x020501A4;
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
      free(aux_types);
      Rprintf("illegal aux_type: %c\n", aux_types[i]);
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
      free( terms[j] );
    }
    free(terms);
  }
  free(aux_types);
  // case then go through and check the values.
  UNPROTECT(1);
  return(ret_data);
}

static const R_CallMethodDef callMethods[] = {
  {"load_bam", (DL_FUNC)&load_bam, 2},
  {"set_iterator", (DL_FUNC)&set_iterator, 3},
  {"clear_iterator", (DL_FUNC)&clear_iterator, 1},
  {"alignments_region", (DL_FUNC)&alignments_region, 8},
  {"sam_read_n", (DL_FUNC)&sam_read_n, 4},
  {"sam_flag_stats", (DL_FUNC)&sam_flag_stats, 2},
  {"target_lengths", (DL_FUNC)&target_lengths, 1},
  {"bam_cigar_str", (DL_FUNC)&bam_cigar_str, 0},
  {"nuc_table", (DL_FUNC)&nuc_table, 0},
  {"bam_flag", (DL_FUNC)&bam_flag, 1},
  {"extract_aux_tags", (DL_FUNC)&extract_aux_tags, 3},
  {NULL, NULL, 0}
};

void R_init_read_bam(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

