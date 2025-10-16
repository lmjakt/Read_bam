#ifndef _QNAME_HASH_H
#define _QNAME_HASH_H

#include "kvec.h"
#include "khash.h"
#include <stdint.h>

// Functions for creating a khash<string, kvec>
// Where kvec holds an array of structs giving the locations of matches

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

typedef struct { kvec_t(sam_record); } sam_record_list;

// A pair of pointers to sam_record_fields held in the hash
// These should probably be considered as volatile if any
// records are added to the sam_record_hash object;
// Hence the data *must* be copied to some other field and the
// user should not modify the hash from which they came.
typedef struct sam_record_pair_p {
  sam_record *read1;
  sam_record *read2;
} sam_record_pair_p;

sam_record_pair_p init_sam_record_pair_p();

// This line will translate into code that defines
// a struct called kh_sam_id_h
// that we should access using the khash_t macro;
// khash_t(sam_id_h) or refer to as sam_rh in klib calls
// 
KHASH_MAP_INIT_STR(sam_id_h, sam_record_list)

// This will translate into a code that defines an integer
// hash. This will be used to keep track of which contigs
// have been indexed. (As we use a map -> vector index
// we need to make sure that we do not index any scaffold
// more than once.
KHASH_SET_INIT_INT(int)

/// A thin wrapper to initialise a sam record hash
/*!    
  @return  A pointer to a khash object of type khash_t(sam_id_h)
*/

khash_t(sam_id_h)* init_sam_hash();

void clear_sam_hash(khash_t(sam_id_h)* hash);

/// Add a sam record to the hash
/*!
  @param rh         A hash of type sam_rh
  @param query_id   The query id
  @param sr         An instance of a sam_record
  @return  An integer giving the number of entries
           in the kvec.
  Note: kh_cstr_t is a typedef of const char*;
*/

int srh_add_element(khash_t(sam_id_h)* hash, kh_cstr_t query_id, sam_record sr);

/// Get a a pointer to a kvec_t(sam_record)
/*!
  @param hash  A pointer to a khash_t object containing query mappings
  @param query_id  The query_id
  @return A pointer to a kvec. 0 if no pointer found.
 */

sam_record_list* srh_get(khash_t(sam_id_h)* hash, kh_cstr_t query_id);

// For hiC type data it is useful to obtain a unique pair of reads and their
// positions. In the simplest case both reads map uniquely to somewhere in the
// assembly. For diploid assemblies, reads can only be uniquely mapping if they
// map to a variant. This complicates the definition of the usable set: we can define
// a scaffold equivalence as a simple integer vector; however, in this case the question
// remains what to return; if we want to return a unique pair, then we have to decide
// which target ids are of interest. We can use the least significant bit of the
// target equivalents to select which of the two will be used; this complicates the API
// but may be worth it to avoid having to pass an additional vector in.

/*! hash: a pointer to 
  @param hash  A pointer to a khash_t object containing query mappings
  @param query_id  The query_id held in a kh_cstr_t object
  @param target_equiv  A pointer to an array giving target_id equivalence mapping; The first
         bit of each value indicates whether reads to the target id represented
	 by the entry (ith position) should be included. The higher bits holds
	 the target id that can be considered as an equivalent. Only one of two
	 equivalents should have the first bit set. If not, and there are alignments
	 to both equivalents then neither will be included.
  @param target_n  The number of entries in target_equiv; this should be equal to the
         total number of targets defined in the bam header.
  @param error: a pointer to an integer; if a unique proper pair is found this will
                be set to 0; otherwise an error code will be set to indicate why
		pair wasn't found.
 */
sam_record_pair_p srh_get_uniq_pair(khash_t(sam_id_h)* hash, kh_cstr_t query_id,
				    uint32_t *target_equiv, size_t target_n,
				    int *error);

#endif
