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

pair_info init_pair_info(sam_record *sr, uint32_t i, uint32_t o, uint32_t read_n, uint32_t mate_n){
  pair_info pi;
  // oo the 'other offset'
  uint32_t oo = (o + 1) % 2;
  uint32_t is_right = (sr[o].target_id > sr[oo].target_id) || (sr[o].begin > sr[oo].begin);
  pi.i = (i << 2) | (o << 1) | is_right;
  pi.pos = sr[o].begin;
  pi.m_target = sr[oo].target_id;
  pi.m_pos = sr[oo].begin;
  pi.al_n = (read_n << 16) | mate_n;
  pi.targets_n = 0;
  return(pi);
}

void set_pair_info(pair_info *pi, sam_record *sr, sam_record *mate_sr, uint32_t i,
		   uint32_t read_n, uint32_t mate_n, uint32_t tg_n1, uint32_t tg_n2){
  //  uint32_t is_right = (sr->target_id > mate_sr->target_id) || (sr->begin > mate_sr->begin);
  // this sets the second bit 
  //uint32_t is_r2 = (sr->flag & 128) >> 6;
  //pi->i = (i << 2) | (is_r2 << 1) | is_right;
  pi->i = i;
  pi->pos = sr->begin;
  pi->m_target = mate_sr->target_id;
  pi->m_pos = mate_sr->begin;
  pi->al_n = (read_n << 16) | mate_n;
  pi->targets_n = (tg_n1 << 16) | tg_n2;
}

target_candidates init_target_candidates(int exempted){
  target_candidates tc;
  kv_init(tc.als[0]);
  kv_init(tc.als[1]);
  tc.min_pos[0] = (uint32_t)-1;
  tc.min_pos[1] = (uint32_t)-1;
  memset(&tc.max_pos, 0, sizeof(tc.max_pos));
  //  memset(&tc.max_score, 0, sizeof(tc.max_score));
  tc.exempted = exempted;
  return(tc);
}

void init_target_candidates_p(target_candidates *tc, int exempted){
  kv_init(tc->als[0]);
  kv_init(tc->als[1]);
  tc->min_pos[0] = (uint32_t)-1;
  tc->min_pos[1] = (uint32_t)-1;
  memset(&tc->max_pos, 0, sizeof(tc->max_pos));
  //  memset(&tc->max_score, 0, sizeof(tc->max_score));
  tc->exempted = exempted;
}

void free_target_candidates(target_candidates *tc){
  kv_destroy(tc->als[0]);
  kv_destroy(tc->als[1]);
}

void clear_target_candidates(target_candidates *tc, size_t i){
  if(i > 1)
    return; // this should be an error.. or warning
  tc->als[i].n = 0;
  tc->min_pos[i] = (uint32_t)-1;
  tc->max_pos[i] = 0;
  //  tc->score[i] = 0;
}

size_t push_target_candidates(target_candidates *tc, size_t i, bam1_t *b){
  if(i > 1)
    return(0);
  sam_record sr;
  sam_record_set(&sr, b);
  kv_push(sam_record, tc->als[i], sr);
  // kv_pushp gives me errors that I do not understand.
  //  sam_record *sr = (sam_record*)kv_pushp( sam_record, tc->als[i] );
  //  sam_record_set( sr, b );
  tc->min_pos[i] = (tc->min_pos[i] > sr.begin) ? sr.begin : tc->min_pos[i];
  tc->max_pos[i] = (tc->max_pos[i] < sr.begin) ? sr.begin : tc->max_pos[i];
  return( tc->als[i].n );
}


//candidate_alignments init_candidate_alignments(size_t n, kvec_t(int) exempted){
candidate_alignments init_candidate_alignments(size_t n, kvi exempted){
  candidate_alignments ca;
  kv_init(ca.targets); // probably not necessary as I then do this manually
  ca.targets.m = n;
  ca.targets.n = n;
  ca.targets.a = malloc(sizeof(target_candidates) * n);
  for(size_t i=0; i < n; ++i)
    init_target_candidates_p(ca.targets.a + i, 0);
  for(size_t i=0; i < exempted.n; ++i){
    if(exempted.a[i] >= 0 && exempted.a[i] < n)
      ca.targets.a[ exempted.a[i] ].exempted = 1;
  }
  for(size_t i=0; i < 4; ++i)
    kv_init(ca.active[i]);
  memset(&ca.score, 0, sizeof(ca.score));
  return(ca);
}

void free_candidate_alignments(candidate_alignments *ca){
  for(size_t i=0; i < ca->targets.n; ++i)
    free_target_candidates(&ca->targets.a[i]);
  kv_destroy(ca->targets);
  for(size_t i=0; i < 4; ++i)
    kv_destroy(ca->active[i]);
}

bit_mask init_bit_mask(size_t n){
  bit_mask bm;
  bm.m = (n % 64 == 0) ? n / 64 : 1 + n/64;
  bm.b_n = n;
  bm.bits = malloc(sizeof(uint64_t) * bm.m);
  memset( bm.bits, 0, sizeof(uint64_t) * bm.m);
  return(bm);
}

// This could be a macro to avoid a function call
void clear_bit_mask(bit_mask *bm){
  memset( bm->bits, 0, sizeof(uint64_t) * bm->m);
}

// This sets bit i to 1 and returns the previous state of i
uint64_t set_bit_mask(bit_mask *bm, size_t i){
  uint64_t b = 1 & (bm->bits[ i / 64 ] >> (i % 64));
  bm->bits[ i / 64 ] |= (1 << (i % 64));
  return(b);
}

void free_bit_mask(bit_mask *bm){
  free(bm->bits);
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
  kv_init(ha.target_n);
  ha.target_bits = init_bit_mask(0);
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
  kv_destroy(ha->target_n);
  free_bit_mask(&ha->target_bits);
}


// read_bits tell us whether we are clearing alignments for read_1
// or read_2.
// read_1 and read_2 have separate vectors for normal and exempted
// scaffolds. The operation (read_bits << 2) | read_bits
// will thus use 4 bits that tell me which bits need to be cleared.
// we also need to clear the max_score for read_1 or read_2
// in order to make a simpler code base. 
void clear_alignments(candidate_alignments *ca, uint32_t read_bits){
  read_bits = (read_bits << 2) | read_bits;
  for(uint32_t k = 0; k < 4; ++k){
    if(read_bits & (1 << k)){
      for(size_t i=0; i < ca->active[k].n; ++i)
	clear_target_candidates( &ca->targets.a[ ca->active[k].a[i] ], k % 2 );
      ca->active[k].n = 0;
      ca->score[ k % 2 ] = 0;
    }
  }
}

void add_alignment(candidate_alignments *ca, bam1_t *b, int al_score){
  size_t read_i = (b->core.flag & 64) ? 0 : 1;
  if(al_score < ca->score[read_i]) // don't add; consider incrementing a counter
    return;
  int tg_id = b->core.tid;
  if(tg_id >= ca->targets.n)
    return;
  if(al_score == ca->score[read_i]){ // add it to the appropriate candidate
    int n = push_target_candidates( ca->targets.a + tg_id, read_i, b);
    if(n == 1)
      kv_push(size_t, ca->active[ read_i + (ca->targets.a[tg_id].exempted * 2) ], tg_id);
    return;
  }
  // A new maximum; clear all the alignments
  clear_alignments(ca, 1 << read_i);
  ca->score[read_i] = al_score;
  push_target_candidates( ca->targets.a + tg_id, read_i, b );
  kv_push(size_t, ca->active[read_i + (ca->targets.a[tg_id].exempted * 2)], tg_id);
}

// NOTE: min_sep defines the minimum distance between read 1 and read 2 of a pair.
//       this separation is not defined if read1 and read2 map to different scaffolds.
//       Given a haplotype resolved assembly, it is to be expected that most read pairs
//       will map to both haplotypes; when any of the potential pairs has a separation
//       that is smaller than min_sep, all of the links should therefore be discarded.
int extract_links(candidate_alignments *ca, uint32_t merge_max_n, uint32_t merge_max_size,
		   uint32_t max_al_n, int min_sep,
		   hiC_assembly *hic, kstring_t query){
  // Allow links between R1 -> R2 (max_al_n each)
  // Allow links between R1 -> X2 (max_al_n, unlimited)
  // Allow links between X1 -> R2
  // Do not allow X1 -> X2;
  // return if alignments to too many non-exempt scaffolds by
  //           read_1 or read_2
  if(ca->active[0].n > max_al_n || ca->active[1].n > max_al_n)
    return(0);
  // return if no alignments to non-exempt scaffolds
  if(ca->active[0].n == 0 && ca->active[1].n == 0)
    return(0);
  // return if alignments by read_1 to both exempt and non-exempt scaffolds
  if(ca->active[0].n > 0 && ca->active[2].n > 0)
    return(0);
  // as above, but for read 2
  if(ca->active[1].n > 0 && ca->active[3].n > 0)
    return(0);
  // return if alignments to exempt scaffolds by read_1 and read_2
  if(ca->active[2].n > 0 && ca->active[3].n > 0)
    return(0);
  // Use exempt or non-exempt scaffolds for read1 and read2
  // ri[0] <- index for read_1, ri[1] <- index for read_2
  size_t ri[2] = {0, 0};
  ri[0] = ca->active[0].n > 0 ? 0 : 2;
  ri[1] = ca->active[1].n > 0 ? 1 : 3;
  //  size_t ri[2] = { ca->active[0].n > 0 ? 0 : 2, ca->active[1].n > 0 ? 1 : 3 };
  // The following check should not really be needed if the checks above do what
  // they should.
  if( ca->active[ri[0]].n == 0 || ca->active[ri[1]].n == 0 )
    return(0);
  // Determine how many will be included
  // This seems wasteful.. 
  //size_t rn[2] = {0,0};   // rn: read number; actually the number of alignments / regions
  //  size_t link_n[2] = {0,0}; // the number of read1 to read2 that should be linked; we don't
                            // know this number beforehand.
                            // lets not bother with determining it first;
                            // lets instead just go through links;
  int count = 0;
  // kv_push is called on:
  //        hic->alignments
  //        hic->query_ids
  //        hic->targets.a[sr_1->target_id].reads
  // we need to be able to revert the push for all of these.
  size_t align_n = hic->alignments.n;
  size_t query_n = hic->query_ids.n;
  // hic->target_n.n is a kvec of size_t pairs (i[2]);
  clear_bit_mask( &hic->target_bits );
  hic->target_n.n = 0;
  int revert = 0;
  // I should probably consider to do this in a different way as revert is likely to
  // end up true for most pairs.
  for(size_t i=0; i < ca->active[ri[0]].n && revert == 0; ++i){
    target_candidates *r1_tg = &ca->targets.a[ ca->active[ri[0]].a[i] ];
    if(r1_tg->exempted == 0 &&
       (r1_tg->als[0].n > merge_max_n || (r1_tg->max_pos[0] - r1_tg->min_pos[0]) > merge_max_size))
      continue;
    // select a single normal read; simply take a median ..
    // custom flag at 4096 to indicate if the record has been recorded.
    int align_i_1 = hic->alignments.n;
    sam_record *sr_1 = &r1_tg->als[0].a[ r1_tg->als[0].n / 2 ];
    for(size_t j=0; j < ca->active[ri[1]].n; ++j){
      target_candidates *r2_tg = &ca->targets.a[ ca->active[ri[1]].a[j] ];
      if(r2_tg->exempted == 0 &&
	 (r2_tg->als[0].n > merge_max_n || (r2_tg->max_pos[0] - r2_tg->min_pos[0]) > merge_max_size))
	continue;
      sam_record *sr_2 = &r2_tg->als[1].a[ r2_tg->als[1].n / 2 ];
      if(abs(sr_1->begin - sr_2->begin) >= min_sep){
	// If min_sep is OK, then create a pair. 
	if((sr_1->flag & 4096) == 0){
	  sr_1->flag |= 4096;
	  //	  align_i_1 = hic->alignments.n;
	  kv_push(sam_record, hic->alignments, *sr_1);
	  kv_push(char*, hic->query_ids, strndup(query.s, query.l));
	}
	pair_info tmp;
	// r1_tg->als[0].n: no. of alignments from the locus that were merged
	// ca->active[ri[0]].n the number of targets to which read 1 mapped
	set_pair_info(&tmp, sr_1, sr_2, align_i_1, r1_tg->als[0].n, r2_tg->als[1].n, ca->active[ri[0]].n, ca->active[ri[1]].n);
	if(!set_bit_mask(&hic->target_bits, sr_1->target_id))
	  kv_push(size_t_pair, hic->target_n, ((size_t_pair){sr_1->target_id, hic->targets.a[sr_1->target_id].reads.n}));
	kv_push(pair_info, hic->targets.a[sr_1->target_id].reads, tmp);
	++count;
	//	kv_push(char*, hic->query_ids, strndup(query.s, query.l));
      }else{
	revert = 1;
	break;
      }
    }
  }
  if(!revert)
    return(count);
  hic->alignments.n = align_n;
  hic->query_ids.n = query_n;
  for(size_t i=0; i < hic->target_n.n; ++i){
    //    Rprintf("setting hic->target_n.a[%ld.i].reads.n (%ld) to %ld\n", hic->target_n.a[i].i,
    //	    hic->targets.a[ hic->target_n.a[i].i ].reads.n,
    //	    hic->target_n.a[i].n);
    hic->targets.a[ hic->target_n.a[i].i ].reads.n = hic->target_n.a[i].n;
  }
  return(0);
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
  hca.target_bits = init_bit_mask(header->n_targets);
  
  size_t count = 0;  // number of pairs extracted
  size_t al_count = 0;  // number of alignments processed
  bam1_t *b = bam_init1();
  if(max_pairs < 1) // consider as unlimited
    max_pairs = (size_t)-1;

  kstring_t last_qid = {0, 0, NULL};
  ksetn_string("blala", 5, &last_qid); // In order to be able to use strcmp;
  //  kstring_t current_qid = {0, 0, NULL};

  // The read_pair holds read one, then read two regardless of the orientation
  // of the alignments.
  sam_record read_pair[2];
  memset(&read_pair, 0, sizeof(sam_record) * 2);

  // use AS (alignment score) as a criterion;
  // we need a max value for read 1 and read 2 of each query.
  int32_t max_AS[2] = {0, 0};
  int32_t AS[2] = {0, 0};
  // Count the number of alignments that have the maximum scores (tie-count)
  uint32_t max_n[2] = {0, 0};
  // I've not found macros for the alignment flag;
  // but the relevant ones are:
  // 4, 8, 64, 256, 512
  // UNMAP, MUNMAP, READ1, SECONDARY DUP
  // WARNING:   we should probably not include 256 in the skip flag here
  //            because the primary alignment can have the same score
  //            as secondary ones.
  uint32_t skip_flag = 4 | 8 | 512;
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
	kv_push( pair_info, hca.targets.a[ read_pair[0].target_id ].reads, init_pair_info(read_pair, hca.alignments.n, 0, max_n[0], max_n[1]) );
	kv_push( sam_record, hca.alignments, read_pair[0] );
	// Then do the same for read 2
	kv_push( pair_info, hca.targets.a[ read_pair[1].target_id ].reads, init_pair_info(read_pair, hca.alignments.n, 1, max_n[1], max_n[0]) );
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
  // free resources.
  bam_destroy1(b);
  sam_hdr_destroy(header);
  sam_close(sam);
  return(hca);
}


// This is the same as extract_read_pairs, except that
// 1. To handle direct repeats that hold very few unique alignments
//    consider up to n alignments separated by a max distance on a single
//    contig as a _single_ alignments. Note that these identical alignments
//    should still have a maximum AS.
// 2. To handle more than one haplotype, allow n regions to hold a maximum
//    alignment score for a given read. In general this will be two. In theory
//    we could have a mapping of contigs considered to be haplotype homologues,
//    but we should not assume such information.
// 3. Allow _special_ exceptions for specified targets. For these exceptions allow
//    all reads that have a maximal score within the specified regions if they
//    have a mate that satisfies the normal criterion.
//    NOTE: the caller must destroy the exemptions vector; it is not done here.
hiC_assembly extract_read_pairs_2(const char* bam_file, int min_sep,
				  int min_AS, size_t max_pairs,
				  int merge_max_n, int merge_max_size,  // -1 for unlimited
				  int max_al_n, kvi exemptions,
				  int *error){
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
  // To simplify this, make a kvec_t, with one entry of the following struct
  // for each target:
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
  hca.target_bits = init_bit_mask(header->n_targets);
  
  size_t count = 0;  // number of pairs extracted
  size_t al_count = 0;  // number of alignments processed
  bam1_t *b = bam_init1();
  if(max_pairs < 1) // consider as unlimited
    max_pairs = (size_t)-1;

  kstring_t last_qid = {0, 0, NULL};
  ksetn_string("blala", 5, &last_qid); // In order to be able to use strcmp;
  //  kstring_t current_qid = {0, 0, NULL};

  candidate_alignments cands = init_candidate_alignments(header->n_targets,
							 exemptions);
  
  // 4, 8, 64, 256, 512
  // UNMAP, MUNMAP, READ1, SECONDARY DUP
  // WARNING:   we should probably not include 256 in the skip flag here
  //            because the primary alignment can have the same score
  //            as secondary ones.
  uint32_t skip_flag = 4 | 8 | 512;
  int align_score = 0;
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

    // if we have a valid pair:
    const char* query_id = bam_get_qname(b);
    if(strncmp(last_qid.s, query_id, last_qid.l)){
      // do important things. This should be a function call as we will need to
      // to do it after the last pair.
      // I should not need to zero the read_pair structure.
      count += extract_links(&cands, (uint32_t)merge_max_n, (uint32_t)merge_max_size, max_al_n, min_sep,
		    &hca, last_qid);
      clear_alignments(&cands, 3);
      ksetn_string(query_id, b->core.l_qname, &last_qid);
    }
    // If the query id is the same as the last one, then continue and
    // look for the maximally scoring alignments;
    // read 1 or read 2
    //size_t i = flag & 64 ? 0 : 1;
    //    int align_score = 0;
    extract_int_aux_values(b, &align_score, (const char*[]){"AS"}, 1);
    if(align_score < min_AS)
      continue;
    add_alignment(&cands, b, align_score);
  }
  count += extract_links(&cands, (uint32_t)merge_max_n, (uint32_t)merge_max_size, max_al_n, min_sep,
			 &hca, last_qid);

  /// Sort all pair_info kvecs; this should be possible to do in separate threads
  /// To speed things up later on.
  for(size_t i=0; i < hca.targets.n; ++i){
    printf("Sorting reads on target %ld\n", i);
    ks_introsort( pair_sort, hca.targets.a[i].reads.n, hca.targets.a[i].reads.a );
  }
  free_candidate_alignments(&cands);
  // free resources associated with sam file operations
  bam_destroy1(b);
  sam_hdr_destroy(header);
  sam_close(sam);
  return(hca);
}
