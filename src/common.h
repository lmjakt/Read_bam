#ifndef _COMMON_H
#define _COMMON_H

#include <R.h>
#include <Rinternals.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <pthread.h>
#include "kvec.h"

// Useful constants:
#define MAX_THREADS 64

// Numeric cigar string operations are defined in sam.h
// mapping 0-9 to:
// BAM_CMATCH, BAM_CINS, BAM_CDEL, BAM_CREF_SKIP, BAM_CSOFT_CLIP, BAM_CHARD_CLIP
// BAM_CPAD, BAM_CEQUAL, BAM_CDIFF, BAM_CBACK 

// cigar operation type defines whether the query, reference
// or both are consumed (1,2,3). 
// sam.h does not define any constants to help us remember which is 1 or 2
// Here, _TP_QC => type consumes query
//       _TP_RC => type consumes reference
#define CIG_OP_TP_QC 1
#define CIG_OP_TP_RC 2

// Bits defining optional fields for:
// alignments_region:

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
#define AR_Q_ERR_EXP 8 // 0x100 estimate the expected number of errors in the sequence.

// A default maximum intron size; this can be over-ridden by the user
// This is relevant to the calculation of intron_depth, but should not
// be used otherwise?
#define MAX_INTRON_L 4096 

// Definitions of data structure fields created by alignments_region
// and sam_read_n

// The number of fields in the list returned by alignments_region()
// and their names.
#define AR_R_FIELDS_N 12
static const char* ar_return_fields[AR_R_FIELDS_N] = {"ref", "query", "al", "ops", "seq", "diff", "depth", "cigar", "qual", "mm", "intron.depth", "mate"};

// The number of rows in an alignments table ($al)
#define AR_AL_RN 15
static const char* ar_al_rownames[AR_AL_RN] = {"flag", "r.beg", "r.end", "q.beg", "q.end", "mqual", "qlen",
					       "qclen", "ops.0", "ops.1", "AS", "NM", "NH", "IH", "exp.err"};
// AR_AL_OPS_ROW specifies the location of ops.0; this is necessary because
// the ops rows need to be set by the cigar parsing function (cigar_to_table)
#define AR_AL_OPS_ROW 8

// the mate quality is specified as an auxiliary tag. We could consider to
// parse the auxiliary information by default but that would be wasteful if
// it is known that it hasn't been set.
#define AR_AL_MATE_RN 3
static const char* al_mate_rownames[AR_AL_MATE_RN] = {"tid", "pos", "isize"};

// The number of rows in an cigar operations table
#define CIG_OPS_RN 8
#define OPS_INIT_SIZE 256
static const char* cig_ops_rownames[CIG_OPS_RN] = {"al.i", "op", "type", "r0", "q0", "r1", "q1", "op.l"};

#define Q_DIFF_RN 4
static const char* q_diff_rownames[Q_DIFF_RN] = {"al.i", "r.pos", "q.pos", "nuc"};

// The fields of base modification information (MM)
#define MM_INFO_RN 7
static const char* mm_info_rownames[MM_INFO_RN] = {"al.i", "q.pos", "mod", "mod.n", "mod.l", "r.pos", "base.inf"};

// The number of rows in a query info table (returned by sam_read_n()) cigar operations table
#define Q_INFO_RN 6
static const char* q_info_rownames[Q_INFO_RN] = {"qlen", "q.cigl", "q.beg", "q.end", "ops.beg", "ops.end"};


// The bit positions used to define what is returned by sam_read_n
// bit 0 indicates the first (least significant bit). To check whether
// bits are set:
// by: 1 << S_QID  => 0, 1 << S_FLAG => 2, etc..
// as implemented in the bit_set() macro (below).
// Note that these are specific to sam_read_n() 
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

// The number of fields returned by sam_read_n
#define SRN_FIELDS_N 16
static const char* srn_return_fields[SRN_FIELDS_N] = {"id", "flag", "ref", "pos", "mapq",
						      "cigar", "ref.m", "pos.m", "tlen", "seq", "qual",
						      "aux", "ops", "q.inf", "mm", "n" };
// The number of standard fields (all vectors of length n):
#define SRN_VEC_FIELDS_N 12
static const unsigned int R_ret_types[SRN_VEC_FIELDS_N] = {STRSXP, INTSXP, INTSXP, INTSXP, INTSXP, // ID, FLAG, RNAME, POS, MAPQ
							   STRSXP, INTSXP, INTSXP, INTSXP, // CIGAR, RNEXT, PNEXT, TLEN
							   STRSXP, STRSXP, STRSXP}; // SEQ, QUAL, AUX

// Nibble -> IUPAC mapping
// Seems like this should be in an htslib header somewhere,
// but I haven't foud it.
// this is a copy of seq_nt16_str defined as a non null-
// terminated vector array
static const char *nuc_encoding = "=ACMGRSVTWYHKDBN";
// if we want to reverse complement
static const char *nuc_encoding_rc = "=TGKCYSBAWRDMHVN";

// The cigar_string struct is probably completely un-necessary
// I suspect that there is an htslib function that I can use to get
// rid of it. But for now leave it here.
struct cigar_string {
  char *buffer;
  size_t string_length; // not including the 0
  size_t buf_size;
};

// Data structures to be used in conjunction with kvec and similar

struct cigar_string cigar_string_init(size_t b_size);

void cigar_string_free( struct cigar_string *cs );

void cigar_string_grow( struct cigar_string *cs );

void cigar_string_read( struct cigar_string *cs, bam1_t *b );

// Structs used by qname_hash and hiC
#define SAM_RECORD_FIELDS_N 10
static const char* sam_record_field_names[SAM_RECORD_FIELDS_N] =
  {"target.id", "r0", "r1", "flag", "q0", "q1", "q.length", "map.q", "qc.length", "AS"};

/*! @typedef 
  @abstract Structure for holding information about a single alignemnt.
  @field target_id    Reference target id
  @field begin        Reference start position
  @field end          Reference end position (optional can be 0)
  @field flag         Sam alignment flag
 */

typedef struct sam_record {
  int target_id;
  int begin;
  int end;
  int flag;
  int q_begin;
  int q_end;
  int q_length;
  int map_q;
  int qc_length;
  int AS; // alignment score; may be NA
} sam_record;

sam_record init_sam_record(int target_id, int begin, int end, int flag,
			   int q_begin, int q_end, int q_length, int map_q, int qc_length,
			   int AS);

void sam_record_set(sam_record *sr, bam1_t *b);

// These are structures for dynamically growing a return data set;
// They would be better off in a library somewhere as I keep repeating
// this code. Or better; use kvec ?

// row major? as in R
struct i_matrix {
  int *data;
  size_t nrow;
  size_t ncol;
  size_t row;
  size_t col;
};

struct i_matrix init_i_matrix(size_t nrow, size_t ncol);
void clear_i_matrix(struct i_matrix *m);

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
struct vectori vectori_init_0(size_t size);
void vectori_push(struct vectori *v, int d);
void vectori_grow_0(struct vectori *v);
void vectori_free(struct vectori *v);
    

// Options taken by cigar_to_table
// depth: if not NULL, pointer for sequence depth
// i_depth; intron_depth
// depth_l: length of depth and i_depth arrays
//
// extra_depth and include_left_als are only relevant to
// multithreaded parsing. See alignments_region_mt_args
// for details.
struct cigar_parse_options {
  int *depth, *i_depth;
  int depth_begin;
  int include_left_als;
  int max_intron_length;
  size_t depth_l;
  char *qseq;
  char *qqual;
  size_t qseq_l;
  const char *rseq;
  size_t rseq_l;
  struct i_matrix *diff, *mm_info;
};

// this defaults to a struct with all elements set to 0;
struct cigar_parse_options init_cigar_parse_options();


// The arguments needed to parse bam files in parallel:
typedef struct alignments_region_mt_args {
  // Variables that should be set by the calling thread:
  int start, end; // the range
  const char *bam_file;
  const char *bam_index_file;
  int target_id;
  struct cigar_parse_options cig_opt;
  size_t imatrix_initial_size;
  uint32_t opts;  // what to return;
  uint32_t *flag_filter; // required and banned flags
  int min_mq; // minimum mapping quality
  int min_ql; // minimum query length;
  ///  size_t extra_depth_isize; // the initial depth size. 

  // Variables that will be initialised in the child thread:
  samFile *sam;
  hts_idx_t *index;
  hts_itr_t *bam_it;   // an initialised bam iterator

  // Data not using the cigar string: these are global and need to be
  // protected by mutex structures when they are updated. Unfortunately
  // They need to be updated in a single go; that means a single mutex
  // to hold them all; Hence prepare all data first and then add to them.
  // If this seems to be a problem, then we can make local copies and do
  // post merging.
  //  pthread_mutex_t *mutex_core;
  struct str_array query_ids;
  struct str_array query_seq;
  struct str_array query_qual;
  struct str_array cigars;
  struct i_matrix al_core; // flag, r.beg, r.edn, q.beg, q.end, mqual, qlen, qclen, ops.0, ops.1
  struct i_matrix mate_core; // tid, pos, isize
  struct i_matrix al_ops;
  // pointers to these will be placed in cig_opt
  struct i_matrix diff, mm_info;
} alignments_region_mt_args;

alignments_region_mt_args init_ar_args();
void free_ar_args(alignments_region_mt_args *args);

typedef struct alignments_merge_args {
  int *al_core;
  int *mate_core;
  int *al_ops;
  int *diff;
  int *mm;

  struct alignments_region_mt_args t_args;
  size_t core_off, ops_off, diff_off, mm_off, mate_off;
} alignments_merge_args;

alignments_merge_args init_ar_merge_args();
void arm_set_offsets(alignments_merge_args *args,
		     size_t core_off, size_t ops_off, size_t diff_off, size_t mm_off, size_t mate_off);

// split_string
// str: the string to be split (0 terminated)
// delim: the char to split by
// n: a pointer to an integer whose value will be set to the number of strings
char **split_string(const char *str, char delim, unsigned int *n);

// Moved from read_bam.c as needed by hiC.c
void extract_int_aux_values(bam1_t *b, int *i_values, const char **tags, size_t tag_n);

// stolen from simple_range/arrange_lines.c
// x1 and x1 are positions of lines that should be drawn
// The function finds y positions that allow this without
// overlaps.
// Note that for each pair of x1 and x2, x1 must be smaller
// than x2
SEXP arrange_lines(SEXP r_x1, SEXP r_x2);

#endif

