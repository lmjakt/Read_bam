#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"

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

sam_record init_sam_record(int target_id, int begin, int end, int flag,
			   int q_begin, int q_end, int q_length, int map_q, int qc_length,
			   int AS){
  sam_record sr;
  sr.target_id=target_id;
  sr.begin=begin;
  sr.end=end;
  sr.flag=flag;
  sr.q_begin=q_begin;
  sr.q_end=q_end;
  sr.q_length=q_length;
  sr.map_q=map_q;
  sr.qc_length = qc_length;
  sr.AS = AS;
  return(sr);
};

void sam_record_set(sam_record *sr, bam1_t *b){
  sr->target_id = (int)b->core.tid;
  sr->begin = (int)b->core.pos;
  sr->flag = (int)b->core.flag;
  sr->q_length = (int)b->core.l_qseq;
  sr->map_q = (int)b->core.qual;
  extract_int_aux_values(b, &sr->AS, (const char*[]){"AS"}, 1);
  // need to parse the cigar values for the others;
  sr->q_begin = 0;
  sr->q_end = 0;
  sr->end = sr->begin;
  sr->qc_length = 0;
  // if no cigar we should not try this
  if(b->core.n_cigar == 0)
    return;
  // parse cigar if present
  uint32_t *cigar = bam_get_cigar(b);
  uint32_t op = bam_cigar_op(*cigar);
  uint32_t op_len = bam_cigar_oplen(*cigar);
  uint32_t op_type = bam_cigar_type(*cigar);
  if(op == BAM_CHARD_CLIP){
    sr->q_begin = op_len;
    sr->q_end = op_len;
    sr->qc_length += op_len;
    if( bam_cigar_op(cigar[1]) == BAM_CSOFT_CLIP )
      sr->q_begin += bam_cigar_oplen(cigar[1]);
  }
  if(op_type & 1){
    sr->q_end = op_len;
    sr->qc_length = op_len;
  }
  if(op_type & 2)
    sr->end += op_len;
  for(int i=1; i < b->core.n_cigar; ++i){
    op_type = bam_cigar_type(cigar[i]);
    op_len = bam_cigar_oplen(cigar[i]);
    sr->q_end += (op_type & 1) ? op_len : 0;
    sr->qc_length += ((op_type & 1) || bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) ?  op_len: 0;
    sr->end += (op_type & 2) ? op_len : 0;
  }
  return;
}

struct i_matrix init_i_matrix(size_t nrow, size_t ncol){
  struct i_matrix m;
  m.nrow=nrow;
  m.ncol=ncol;
  m.row=0;
  m.col=0;
  m.data = malloc( sizeof(int) * m.nrow * m.ncol );
  return(m);
}

void clear_i_matrix(struct i_matrix *m){
  m->nrow=0;
  m->ncol=0;
  m->row=0;
  m->col=0;
  free(m->data);
  m->data = 0;
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

// copies n bytes; terminates with a 0
void str_array_push_cp_n(struct str_array *str, const char *word, size_t n){
  if(str->length >= str->capacity)
    double_str_array(str);
  str->strings[str->length] = malloc(sizeof(char) * (n + 1));
  str->strings[str->length][ n ] = 0;
  memcpy( str->strings[str->length], word, sizeof(char) * (n));
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


struct vector vector_init(size_t capacity, size_t unit_size){
  struct vector v;
  v.capacity = capacity;
  v.length = 0;
  v.unit_size = unit_size;
  v.data = malloc( capacity * unit_size );
  return(v);
}

void vector_push(struct vector *v, void *data){
  if(v->length + 1 > v->capacity){
    v->capacity *= 2;
    v->data = realloc(v->data, v->capacity * v->unit_size);
  }
  memcpy( v->data + v->length, data, v->unit_size );
}

void* vector_at(struct vector *v, size_t i){
  if(i < v->length)
    return( v->data + i * v->unit_size );
  return(0);
}

void vector_free(struct vector *v){
  free(v->data);
  v->data = 0;
  v->capacity = 0;
  v->length = 0;
}

void vector_clear(struct vector *v){
  v->length = 0;
}


struct vectori vectori_init(size_t size){
  struct vectori v;
  v.n = 0;
  v.m = size;
  v.data = (v.m > 0) ? malloc(v.m * sizeof(int)) : 0;
  return(v);
}

void vectori_push(struct vectori *v, int d){
  if(v->n >= v->m){
    v->m = (v->m == 0) ? 1 : v->m << 1;
    v->data = realloc(v->data, sizeof(int) * v->m);
  }
  v->data[v->n] = d;
  v->n++;
}

void vectori_grow_0(struct vectori *v){
  size_t om = v->m;
  v->m = 2 * v->m;
  v->data = realloc(v->data, sizeof(int) * v->m);
  memset(v->data + om, 0, sizeof(int) * om);
}

void vectori_free(struct vectori *v){
  free(v->data);
  v->data = 0;
  v->m = 0;
  v->n = 0;
}


// this defaults to a struct with all elements set to 0;
struct cigar_parse_options init_cigar_parse_options(){
  struct cigar_parse_options cpo;
  memset( &cpo, 0, sizeof(struct cigar_parse_options) );
  // include_left_als is usually allowed
  cpo.include_left_als = 1;
  return(cpo);
}

alignments_region_mt_args init_ar_args(){
  alignments_region_mt_args args;
  memset(&args, 0, sizeof(alignments_region_mt_args));
  args.cig_opt = init_cigar_parse_options();
  args.imatrix_initial_size = OPS_INIT_SIZE;
  return(args);
}

alignments_merge_args init_ar_merge_args(){
  alignments_merge_args args;
  memset(&args, 0, sizeof(alignments_merge_args));
  return(args);
}

void free_ar_args(alignments_region_mt_args *args){
  str_array_free( &args->query_ids );
  str_array_free( &args->query_seq );
  str_array_free( &args->query_qual );
  str_array_free( &args->cigars );
  clear_i_matrix( &args->al_core );
  clear_i_matrix( &args->mate_core );
  clear_i_matrix( &args->al_ops );
  clear_i_matrix( &args->diff );
  clear_i_matrix( &args->mm_info );
}

void arm_set_offsets(alignments_merge_args *args,
		     size_t core_off, size_t ops_off, size_t diff_off, size_t mm_off, size_t mate_off){
  args->core_off = core_off;
  args->ops_off = ops_off;
  args->diff_off = diff_off;
  args->mm_off = mm_off;
  args->mate_off = mate_off;
}

// This allocates new words; that may or may not be what you want
// to do;
void set_word(char **words, const char* beg, const char* end){
  *words = malloc(1 + end - beg);
  strncpy( *words, beg, end-beg );
  (*words)[end-beg] = 0;
}

char** split_string(const char *str, char delim, unsigned int *n){
  unsigned int str_l = 0;
  *n = 1;
  //  unsigned int word_count = 1;
  const char *end = str;
  const char *beg = str;
  while(*end){
    if(*end == delim)
      (*n)++;
    str_l++;
    ++end;
  }
  char **words = malloc(sizeof(const char*) * (*n));
  end = str;
  unsigned int word_i = 0;
  while(*end){
    if(*end == delim){
      set_word( words + word_i, beg, end );
      beg = end+1;
      ++word_i;
    }
    ++end;
  }
  set_word( words + word_i, beg, end );
  return(words);
}

// Moved from read_bam.c to be available to other compilation units.
void extract_int_aux_values(bam1_t *b, int *i_values, const char **tags, size_t tag_n){
  for(size_t i=0; i < tag_n; ++i){
    i_values[i] = R_NaInt;
    uint8_t *s = bam_aux_get(b, tags[i]);
    if(!s)
      continue;
    errno = 0;
    int64_t v = bam_aux2i(s);
    i_values[i] = (errno != EINVAL) ? (int)v : R_NaInt;
  }
}


// The following function was copied from R/simple_range/arrang_lines.c
// It doesn't really belong here, but it's useful for visualising alignments
// I may want to move this to an util.c or similar file instead. But for now
// lets stick it here and see if I can get it to work.

// x and y are positions of lines that need to be drawn
// The function finds y positions that allow this
// Note that for each pair of x1 and x2, x1 must be smaller
// than x1
SEXP arrange_lines(SEXP r_x1, SEXP r_x2){
  if(!isReal(r_x1) || !isReal(r_x2))
    error("Both arguments should be real values\n");
  if(length(r_x1) != length(r_x2) || length(r_x1) == 0)
    error("Arguments should have the same non-zero length\n");
  int n = length(r_x1);
  double *x1 = REAL(r_x1);
  double *x2 = REAL(r_x2);
  // assign a vector of ints..
  SEXP r_y = PROTECT(allocVector(INTSXP, n));
  int *y = INTEGER(r_y);
  // use a char array as a boolean vector. It's ugly, but
  // it seems to work. We probably have a bit of a slow
  // down due to the repeated calls to bzero. There should
  // be a better way of handling it... We could easily use
  // bit bucket as well, but that's probably more trouble than it's
  // worth.
  unsigned char *forbidden_y = malloc(n); // use as a boolean vector 
  bzero((void*)y, sizeof(int) * n);
  bzero((void*)forbidden_y, n);
  for(int i=1; i < n; ++i){
    int max_y = 0;
    for(int j=0; j < i; ++j){
      if( x2[j] >= x1[i] && x2[i] >= x1[j]){
	forbidden_y[y[j]] = 1;
	max_y = y[j] > max_y ? y[j] : max_y;
      }
    }
    for(int j=0; j < i; ++j){
      if(!forbidden_y[j]){
	y[i] = j;
	break;
      }
      // if no empty slot found, set to position to i
      y[i] = (y[i] == 0) ? i : y[i];
    }
    bzero((void*)forbidden_y, n);
  }
  free(forbidden_y);
  UNPROTECT(1);
  return(r_y);
}

