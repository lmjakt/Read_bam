#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"


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

void vectori_free(struct vectori *v){
  free(v->data);
  v->data = 0;
  v->m = 0;
  v->n = 0;
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

