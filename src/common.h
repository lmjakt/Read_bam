#ifndef _COMMON_H
#define _COMMON_H

#include <R.h>
#include <Rinternals.h>

// Useful constants:

// The following are true for the encoding used in cig_ops tables;
// as provided to R from functions here.
// They differ from the original bam encoding by counting from 0.
// it's BAM encoding + 1
// We probably don't need these; we should instead use the
// BAM_CMATCH, BAM_CINS, BAM_CDEL, BAM_CREF_SKIP, BAM_CSOFT_CLIP, BAM_CHARD_CLIP
// BAM_CPAD, BAM_CEQUAL, BAM_CDIFF, BAM_CBACK that are defined in sam.h
#define CIG_M 1
#define CIG_I 2
#define CIG_D 3
#define CIG_N 4
#define CIG_S 5
#define CIG_H 6
#define CIG_P 7
#define CIG_EQ 8
#define CIG_X 9

// cigar operation type defines whether the query, reference
// or both are consumed (1,2,3). 
// sam.h does not define any constants to help us remember which is 1 or 2
// Here, _TP_QC => type consumes query
//       _TP_RC => type consumes reference
#define CIG_OP_TP_QC 1
#define CIG_OP_TP_RC 2

// The following define the number of rows  and their names in tables returned by
// aligned_region and sam_read_n

// The number of rows in an alignments table ($al)
#define AR_AL_RN 10
static const char* ar_al_rownames[AR_AL_RN] = {"flag", "r.beg", "r.end", "q.beg", "q.end", "mqual", "qlen", "qclen", "ops.0", "ops.1"};

// The number of rows in an cigar operations table
#define CIG_OPS_RN 8
#define OPS_INIT_SIZE 256
// a global variable; I feel dirty;
static const char* cig_ops_rownames[CIG_OPS_RN] = {"al.i", "op", "type", "r0", "q0", "r1", "q1", "op.l"};

#define Q_INFO_RN 6
static const char* q_info_rownames[Q_INFO_RN] = {"qlen", "q.cigl", "q.beg", "q.end", "ops.beg", "ops.end"};

#define MM_INFO_RN 7
static const char* mm_info_rownames[MM_INFO_RN] = {"al.i", "q.pos", "mod", "mod.n", "mod.l", "r.pos", "base.inf"};


// not sure why this isn't in a header somewhere
// but I can't find it.
// this is a copy of seq_nt16_str defined as a non null-
// terminated vector array
static const char *nuc_encoding = "=ACMGRSVTWYHKDBN";
// if we want to reverse complement
static const char *nuc_encoding_rc = "=TGKCYSBAWRDMHVN";

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

struct i_matrix init_i_matrix(size_t nrow, size_t ncol);

void double_columns( struct i_matrix *m );


// push numbers in c onto the matrix
void push_column( struct i_matrix *m, int *c );

struct str_array {
  size_t capacity;
  size_t length;
  char **strings;
};

struct str_array init_str_array(size_t init_size);

void double_str_array(struct str_array *str);

// copy the 0 delimited word to str_array
void str_array_push_cp(struct str_array *str, const char *word);
// copies n bytes; terminates with a 0
void str_array_push_cp_n(struct str_array *str, const char *word, size_t n);
// add pointer to word without copying.
void str_array_push(struct str_array *str, char *word);
void str_array_free(struct str_array *str);

// I suspect that libc may have some more efficient
// implementation of this:
struct vector {
  size_t capacity;
  size_t length;
  size_t unit_size;
  void *data;
};

struct vector vector_init(size_t capacity, size_t unit_size);

void vector_push(struct vector *v, void *data);
void* vector_at(struct vector *v, size_t i);
void vector_free(struct vector *v);
void vector_clear(struct vector *v);

// More reimplementation;
// This should be moved to klib; kvec etc which
// provide more consistent and convenient interfaces.

struct vectori {
  int *data;
  size_t n; // the number of elements stored
  size_t m; // the capacity
};

// Functions to initialise push and free;
struct vectori vectori_init(size_t size);
void vectori_push(struct vectori *v, int d);
void vectori_free(struct vectori *v);
    


// split_strong
// str: the string to be split (0 terminated)
// delim: the char to split by
// n: a pointer to an integer whose value will be set to the number of strings
char **split_string(const char *str, char delim, unsigned int *n);

// stolen from simple_range/arrange_lines.c
// x1 and x1 are positions of lines that should be drawn
// The function finds y positions that allow this without
// overlaps.
// Note that for each pair of x1 and x2, x1 must be smaller
// than x2
SEXP arrange_lines(SEXP r_x1, SEXP r_x2);

#endif

