#include "kstring.h"
#include "ksort.h"
#include <htslib/sam.h>
#include <htslib/hts.h>
#include "common.h"
#include "hiC.h"

// Note use of ksort is as described at:
// http://attractivechaos.github.io/klib/#Ksort%3A%20sorting%2C%20shuffling%2C%20heap%20and%20k-small
//
// Discussion from 2008 / 2009 at:
// https://attractivechaos.wordpress.com/2008/08/28/comparison-of-internal-sorting-algorithms/
//
// Recommends introsort when comparison is cheap; also suggests that glibc qsort is
// slow. Whether that is still true is questionable.

// A macro to define an array of pair_info structs
// defined in hiC.h
#define pair_info_lt(a, b) ((a).pos < (b).pos)

KSORT_INIT(pair_sort, pair_info, pair_info_lt);

pair_info init_pair_info(sam_record *sr, uint32_t i, uint32_t o){
  pair_info pi;
  // oo the 'other offset'
  uint32_t oo = (o + 1) % 2;
  uint32_t is_right = (sr[o].target_id > sr[oo].target_id) || (sr[o].begin > sr[oo].begin);
  pi.i = (i << 2) | (o << 1) | is_right;
  pi.pos = sr[o].begin;
  pi.m_target = sr[oo].target_id;
  pi.m_pos = sr[oo].begin;
  return(pi);
}

void init_hiC_target(hiC_target *ht, int32_t target_id, int32_t target_length){
  ht->target_id = target_id;
  ht->target_length = target_length;
  kv_init(ht->reads);
}

void free_hiC_target(hiC_target *ht){
  kv_destroy(ht->reads);
}

hiC_assembly init_hiC_assembly(){
  hiC_assembly ha;
  kv_init(ha.target_lengths);
  kv_init(ha.targets);
  kv_init(ha.alignments);
  kv_init(ha.query_ids);
  return(ha);
}

void free_hiC_assembly(hiC_assembly *ha){
  kv_destroy(ha->target_lengths);
  for(size_t i=0; i < ha->targets.n; ++i)
    free_hiC_target(ha->targets.a + i);
  kv_destroy(ha->targets);
  for(size_t i=0; i < ha->query_ids.n; ++i)
    free( ha->query_ids.a[i] );
  kv_destroy(ha->query_ids);
  kv_destroy(ha->alignments);
}

hiC_assembly extract_read_pairs(const char* bam_file, int min_sep, int min_qual,
				int min_AS, size_t max_pairs, int *error){
  *error = 1;
  hiC_assembly hca = init_hiC_assembly();
  samFile *sam = sam_open( bam_file, "r" );
  if(!sam){
    *error = 1;
    return(hca);
  }
  sam_hdr_t *header = sam_hdr_read(sam);
  if(!header){
    *error = 2;
    sam_close(sam);
    return(hca);
  }
  // we can copy the data from the sam_hdr_t to the target lengths vector;
  // we keep this in a kvec (though it's not really necessary as we define the
  // number of entries at this stage).
  hca.target_lengths.m = header->n_targets;
  hca.target_lengths.n = header->n_targets;
  hca.target_lengths.a = realloc( hca.target_lengths.a, sizeof(uint32_t) * header->n_targets );
  memcpy( hca.target_lengths.a, header->target_len, sizeof(uint32_t) * header->n_targets );
  // similarly we can initialise the new data for the targets
  hca.targets.m = header->n_targets;
  hca.targets.a = realloc( hca.targets.a, sizeof(hiC_target) * header->n_targets );
  hca.targets.n = header->n_targets;
  for(int32_t i=0; i < header->n_targets; ++i)
    init_hiC_target(hca.targets.a + i, i, header->target_len[i]);
  
  size_t count = 0;  // number of pairs extracted
  size_t al_count = 0;  // number of alignments processed
  bam1_t *b = bam_init1();
  if(max_pairs < 1) // consider as unlimited
    max_pairs = (size_t)-1;

  kstring_t last_qid = {0, 0, NULL};
  ksetn_string("blala", 5, &last_qid); // In order to be able to use strcmp;
  kstring_t current_qid = {0, 0, NULL};

  // The read_pair holds read one, then read two regardless of the orientation
  // of the alignments.
  sam_record read_pair[2];
  memset(&read_pair, 0, sizeof(sam_record) * 2);

  // use AS (alignment score) as a criterion;
  // we need a max value for read 1 and read 2 of each query.
  uint32_t max_AS[2] = {0, 0};
  uint32_t AS[2] = {0, 0};
  // Count the number of alignments that have the maximum scores (tie-count)
  uint32_t max_n[2] = {0, 0};
  // I've not found macros for the alignment flag;
  // but the relevant ones are:
  // 4, 8, 64, 256, 512
  // UNMAP, MUNMAP, READ1, SECONDARY DUP
  uint32_t skip_flag = 4 | 8 | 256 | 512;
  while( count < max_pairs && sam_read1(sam, header, b) >= 0){
    ++al_count;
    if(al_count % 10000 == 0){
      printf("Processed %ld alignments to %ld pairs\n", al_count, count);
    }
    uint32_t flag = b->core.flag;
    if(flag & skip_flag)
      continue;
    // get query id and check if query id is different from previous
    // query id. 

    // If the query id is not the same as the last one, then check
    // if we have a valid pair:
    const char* query_id = bam_get_qname(b);
    if(strncmp(last_qid.s, query_id, last_qid.l)){
      // do important things. This should be a function call as we will need to
      // to do it after the last pair.
      if(max_n[0] == 1 && max_n[1] == 1 && (read_pair[0].target_id != read_pair[1].target_id ||
					   abs(read_pair[0].begin - read_pair[1].begin) >= min_sep)){
	// Push read 1 info to the target struct and then the full sam record
	kv_push( pair_info, hca.targets.a[ read_pair[0].target_id ].reads, init_pair_info(read_pair, hca.alignments.n, 0) );
	kv_push( sam_record, hca.alignments, read_pair[0] );
	// Then do the same for read 2
	kv_push( pair_info, hca.targets.a[ read_pair[1].target_id ].reads, init_pair_info(read_pair, hca.alignments.n, 1) );
	kv_push( sam_record, hca.alignments, read_pair[1] );
	kv_push( char*, hca.query_ids, strndup(last_qid.s, last_qid.l) );
	++count;
      }
      // we can't unfortunately use AS = {0,0}; that's only ok for initialisation
      memset(AS, 0, sizeof(AS));
      memset(max_AS, 0, sizeof(max_AS));
      memset(max_n, 0, sizeof(max_n));
      ksetn_string(query_id, b->core.l_qname, &last_qid);
      // I should not need to zero the read_pair structure. 
    }
    // If the query id is the same as the last one, then continue and
    // look for the maximally scoring alignments;
    // read 1 or read 2
    size_t i = flag & 64 ? 0 : 1;
    extract_int_aux_values(b, (AS + i), (const char*[]){"AS"}, 1);
    if(AS[i] == max_AS[i])
      max_n[i]++;
    if(AS[i] > max_AS[i] && AS[i] >= min_AS && b->core.qual >= min_qual){
      max_AS[i] = AS[i];
      max_n[i] = 1;
      // and obtain the read information.
      sam_record_set( read_pair + i, b );
    }
  }
  // At this point we should add the last read pair to the data; however, I will
  // leave that until I have implemented a function to do this as it will make the
  // code horrible if I do that now. First check that the code above works and then
  // refactor it to fix this issue.

  /// Sort all pair_info kvecs; this should be possible to do in separate threads
  /// To speed things up later on.
  for(size_t i=0; i < hca.targets.n; ++i){
    printf("Sorting reads on target %ld\n", i);
    ks_introsort( pair_sort, hca.targets.a[i].reads.n, hca.targets.a[i].reads.a );
  }
  return(hca);
}
