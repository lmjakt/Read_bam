#ifndef _QNAME_HASH_H
#define _QNAME_HASH_H

#include "kvec.h"
#include "khash.h"

// Functions for creating a khash<string, kvec>
// Where kvec holds an array of structs giving the locations of matches

#define SAM_RECORD_FIELDS_N 9
static const char* sam_record_field_names[SAM_RECORD_FIELDS_N] =
  {"target.id", "r0", "r1", "flag", "q0", "q1", "q.length", "map.q", "qc.length"};

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
} sam_record;

sam_record init_sam_record(int target_id, int begin, int end, int flag,
			   int q_begin, int q_end, int q_length, int map_q, int qc_length);

typedef struct { kvec_t(sam_record); } sam_record_list;

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
  @param query_id  The query_id
  @return A pointer to a kvec. 0 if no pointer found.
 */

sam_record_list* srh_get(khash_t(sam_id_h)* hash, kh_cstr_t query_id);

#endif
