#include "kvec.h"
#include "khash.h"
#include "qname_hash.h"

sam_record init_sam_record(int target_id, int begin, int end, int flag,
			   int q_begin, int q_end, int q_length, int map_q, int qc_length){
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
  return(sr);
};


khash_t(sam_id_h)* init_sam_hash(){
  return(kh_init(sam_id_h));
};


// this needs to iterate over all entries in order to clear
// the kvec_t objects
void clear_sam_hash(khash_t(sam_id_h)* hash){
  khiter_t k;
  for(k = kh_begin(hash); k != kh_end(hash); ++k){
    if(kh_exist(hash, k)){
	kv_destroy( kh_val(hash, k) );
	kh_del( sam_id_h, hash, k );
    }
  } 
}


// empty stub for now to see if I'm getting the declarations correct..
int srh_add_element(khash_t(sam_id_h) *hash, kh_cstr_t query_id, sam_record sr){
  // first we use kh_get; this gives an iterator if an entry exists
  // we then use kh_value, and append or make a new kvec_t and
  // put that into the hash.
  int ret = 0;
  khiter_t k = kh_get(sam_id_h, hash, query_id);
  if(k == kh_end(hash)){
    k = kh_put(sam_id_h, hash, strdup(query_id), &ret);
    // might be an idea to check the value of k here
    if(k == kh_end(hash))
      return(0);
    kv_init(kh_val(hash, k));
  }
  kv_push( sam_record, kh_val(hash, k), sr );
  return(1);
}

sam_record_list* srh_get(khash_t(sam_id_h) *hash, kh_cstr_t query_id){
  khiter_t k = kh_get(sam_id_h, hash, query_id);
  if( k == kh_end(hash) )
    return(0);
  return(&kh_val(hash, k));
}


