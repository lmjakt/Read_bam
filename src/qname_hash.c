#include "kvec.h"
#include "khash.h"
#include "qname_hash.h"
#include <stdio.h>



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


sam_record_pair_p init_sam_record_pair_p(){
  sam_record_pair_p srp;
  srp.read1 = 0;
  srp.read2 = 0;
  return(srp);
}

sam_record_pair_p srh_get_uniq_pair(khash_t(sam_id_h)* hash, kh_cstr_t query_id,
				    uint32_t *target_equiv, size_t target_n, int *error){
  sam_record_list* srl= srh_get(hash, query_id);
  sam_record_pair_p srp = init_sam_record_pair_p();
  // Since only unique pairs are allowed, and since each read can only map
  // to two different locations, the maximum allowed number of alignments
  // is 4. Any more than that and we should simply return an empty set;
  // Similarly we should have at least two alignments;
  if(srl->n > 4 || srl->n < 2){
    *error = 1;
    return(srp);
  }

  // the sam record list indices for read1 and read2
  size_t read1_i[2] = {0,0};
  size_t read2_i[2] = {0,0};
  // and the taget ids to which the reads have been aligned
  size_t read1_target[2] = {0, 0};
  size_t read2_target[2] = {0, 0};

  // The number of alignments for read 1 and 2. Both should
  // be 1 or 2;
  int read1_n = 0;
  int read2_n = 0;

  // bitwise masks holding which of the two possible reads
  // are allowed; only one for each should be;
  unsigned int r1_b, r2_b;
  r1_b = r2_b = 0;

  // then go through the entries and see if we can define both read1 and read2
  // appropriately:
  /// This needs to be modified such that map quality is taken into consideration
  /// It would be better to use something that allows for two top scoring alignments
  /// but at the moment that is likely to be difficult.
  for(size_t i=0; i < srl->n; ++i){
    if(srl->a[i].target_id >= target_n){ /// this should really warn the user
      *error = 2;
      return(srp);
    }
    int is_read1 = srl->a[i].flag & 64;
    if( (is_read1 && read1_n > 1) || (!is_read1 && read2_n > 1)){
      *error = 3;
      return(srp);
    }
    if(is_read1){
      read1_i[ read1_n ] = i;
      read1_target[ read1_n ] = srl->a[i].target_id;
      r1_b |= ((target_equiv[read1_target[read1_n]] & 1) ? 1 << read1_n : 0);
      read1_n++;
    }
    if(!is_read1){
      read2_i[ read2_n ] = i;
      read2_target[ read2_n ] = srl->a[i].target_id;
      r2_b |= ((target_equiv[read2_target[read2_n]] & 1) ? 1 << read2_n : 0);
      read2_n++;
    }
  }
  // Check that we have exactly one appropriate target for each
  // read:
  // the number of set bits should be 1 for both; i.e. values of 1 or 2
  if( __builtin_popcount(r1_b) != 1 || __builtin_popcount(r2_b) != 1 ){
    *error = 4;
    return(srp);
  }

  // if we are here then it is established that there is exactly
  // one acceptable alignment for read1 and read2; We only need
  // to confirm that the alignments are to equivalent targets if
  // there is more than one:
  // THE equivalents vector must be using a 1 based offset; 0 should indicate
  // no homologous sequence; to compensate add one to the target obtained.
  if(read1_n > 1 && (target_equiv[ read1_target[0] ] >> 1) != (1+read1_target[1])){
    *error = 5;
    return(srp);
  }
  if(read2_n > 1 && (target_equiv[ read2_target[0] ] >> 1) != (1+read2_target[1])){
    *error = 6;
    return(srp);
  }

  *error = 0;
  srp.read1 = &(srl->a[ read1_i[ (r1_b == 1) ? 0 : 1 ] ]);
  srp.read2 = &(srl->a[ read2_i[ (r2_b == 1) ? 0 : 1 ] ]);
  return(srp);
}
