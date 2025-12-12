#ifndef _HIC_H
#define _HIC_H

#include <stdint.h>
#include "kstring.h"
#include "kvec.h"
#include "common.h"

////// These functions are meant to consider how to:
////// 1. Compensate for uneven coverage of uniqualy mapping reads in HiC-data
////// 2. Infer distances and (thus) orientations between contigs / scaffolds
///
/// The equations used are based on the observation that cumulative number
/// of (distant) read-pairs within a distance d increases linearly with log(d)
/// This implies (?) that that p[d] ~ 1/d
/// since the integral of 1/d is log(d)
///
/// This means that if we have two regions A and B of length l_a and l_b
/// and where d is defined as l_a + d_sep (the distance separating the two
/// regions) and n1 and n2 are denoted by the number of reads where
/// both read1 and read1 map to A (n1), or where the reads link A and B (n2).
///
/// n2/n1 = f = ((log(l_b + d) - log(d)) / (log(l_a))) * (N_b / N_a)
///
/// Where N_a and N_b are the total number of reads that map to regions A and
/// B. This should compensate for mapping density and random links.
///
/// This equation can be rearranged to give d from n2/n1 (f)
///
/// set K (a known quantity) : 
/// K <- (f * log(l_a)) / (N_b / N_a)
/// then:
///
/// log(l_b + d) - log(d) = K
///
/// log((l_b + d) / d) = K
///
/// (l_b + d) / d = exp(K)
///
/// l_b + d = d * exp(K)
///
/// l_b + d - d * exp(K) = 0
///
/// l_b + d(1 - exp(K)) = 0
///
/// d(1 - exp(K)) = -l_b
///
/// d = -l_b / (1 - exp(K))
///
/// d = l_b / (1 - exp(K))
///
/// If we expand K:
///
/// d = l_b / (1 - exp( ((n2/n1) * log(l_a)) / (N_b / N_a) ))
///
/// This can trivially be tested using scaffolds that are thought to
/// be correct.
/// Note that the equation should possibly also be adjusted so that it compensates
/// for the fact that I don't allow read pairs with isizes below some limit
/// I'm not quite sure right now of the best way to do this though.

// Approximate procedure to be followed:
// 1. Read in bam file ordered by read-name;
//    a. Identify uniquely mapped pairs having suitable properties (eg.
//       uniquely mapping and with sufficient separation).
//    b. Assign to target structs with information about mates
// 2. Sort alignments on each target independently
// 3. For each target define contiguous regions with similar map densities
//    and record these in the structs.
// 4. Use these regions to estimate interegion distances measures.
// 5. Evaluate result of distances; consider visualisation and rearrangement
//    if evidence for additional linkage (re-scaffolding). How to do this
//    I'm not sure, but there are many approaches.


// Functions for obtaining read-pairs that are suitable for HiC analyses
//
// There are two options:
// Option 1:
// The information about read pair mappings is split into a set of structs
// Each struct contains only integer, double or char* data;
// Data describing read pairs will be stored as kvec vectors of such structs.
// Restricting them to hold only a single data type means that these can
// be directly converted to R matrix or vector objects; either by:
// memcpy() to R data
// or
// by using allocVector3 with a custom allocator function.
//
// Option 2:
// Keep all information in a kvec vector of a single struct type; the data elements
// of this struct will then have to be copied to R data structures by looping through
// each element and doing large numbers of assignment operations.
//
// Option 3:
// Define a struct for a single alignment of the pair;
//        Keep an array of two members [2];
//        two integers which define which is the left and which is the right
//        and a char* for the query id;
//
// Option 4:
// All read pairs are kept in a kvec(hiC_read_pair) and are ordered low to high
//     based on the ordering of target sequences and positions within them;
//     i.e. the information is kept as a:
//     [[left], [right]], [[left], [right]], ...
//     This is a global vector containing all alignment information in the file or files.
//
// There is also a corresponding global vector of kstring_ts
// kvec(kstring_t) with each entry corresponding to a read pair
//
// Information about the alignments to a target sequence is then kept as:
// kvec(size_t) left
// kvec(size_t) right
// and this can be reordered if needed. We can also keep a logical value indicating
// whether the order of alignments along the target should be considered inverted
// and whether (separately) the pairs should be considered inverted. If that is the case
// then left and right reads should be considered separately.
// 
// 
// 
// Option 1 requires more complexity in the phase of finding the read pairs, and will
// lead to a larger number of memory allocations during the data parsing operation.
// This is likely to be the bottleneck of the process. 
//
//
// Option 2 is simpler during the parsing phase, but is more complex in the creation
// of R data structures. In particular it will make it more cumbersome to modify the
// datastructures in the future.
//
// Option 3 would simplify the swapping of left and right positions if we want to reverse
//        the orientation of targets.
// 
// 
// Option 4 allows us to keep links of both right and left reads from a given chromosome.
//          without using up a lot more data (64 bits for each). This is thus the simplest
//          option. The sam_record is also a simple, all int data structure and so it may
//          be possible to transfer to R with a single memcpy command.
//
// I'm now using a variant of option 4. However, In addition to the index in the global
// vector of sam_records the target vectors also keep some minimal information about position
// and mate position (including other target). This should make for faster traversall of data
// but breaks the golden rule of data structures by keeping some information in three places at
// the same time. This arrangement of data should thus be considered as experimental.




// sam_record is defined in common.h
// 
// all the fields in sam_record are simple ints; this may cause problems
// in the future, and it also means that we need to cast int64_t to int32_t
// as the bam structure uses 64 bits for the position.
// This is not ideal, but for now lets keep it and monitor memory consumption.


/*! @typedef
  @abstract Holds minimal data about a pair of alignments; only the mate_target is specified
            as the pair_info struct will be held in a vector bound to a specific target.
  @field    i         upper 30 bits hold the index to a global kvec of sam_records
		      bit 2: 0 if read is read one
		      bit 1: 0 if read is left read
  @field    pos       position of read (reference begin by default)
  @field    m_target  target_id of mate mapping
  @field    m_pos     position of mate
  @field    al_n      The number of alignments of read and mate with the same alignment score
                      held in the lower and upper 16 bits for the mate and read respectively

  @Discussion  Using only 30 bits to hold the index means that we are limited to using
               about 500 million read pairs. That *should* generally be OK. But it should not
	       be difficult to change this if needed. The struct is in any case experimental
	       and I expect it will change with time. I'm not really sure whether we need to
	       define left and right in the struct itself.
 */
typedef struct {
  uint32_t i;
  uint32_t pos;
  uint32_t m_target;
  uint32_t m_pos;
  uint32_t al_n;
  uint32_t targets_n;
} pair_info;


// Macros to access elements of pair_info::i
#define pi_is_left(p) (((p.i) & 1) == 0)
#define pi_is_r1(p) (((p.i) & 2) == 0)
#define pi_i(p) ( (p.i) >> 2 )

// sr should be sam_record sr[2]; not sure I can define that
// in the header.
pair_info init_pair_info(sam_record *sr, uint32_t i, uint32_t o, uint32_t read_n, uint32_t mate_n);

// Holds alignments for members of a read pair that have a maximum score
// and are on the same target
//typedef kvec_t(sam_record) sam_record_list;

typedef kvec_t(int) kvi;

typedef struct {
  kvec_t(sam_record) als[2];
  uint32_t min_pos[2];
  uint32_t max_pos[2];
  //  int score[2];
  int exempted;  // 0 for normal; 1 for exemptions.
} target_candidates;

target_candidates init_target_candidates(int exempted);
void init_target_candidates_p(target_candidates *tc, int exempted);
void free_target_candidates(target_candidates *tc);

typedef struct {
  kvec_t(target_candidates) targets;
  // candidates with max score on read 1 or read 2, normal or exempt
  // norm read1, norm read2, exempt read1, exempt read2
  kvec_t(size_t) active[4]; 
  int score[2]; // Alignment score for read 1 and read2
} candidate_alignments;

//candidate_alignments init_candidate_alignments(size_t n, kvec_t(int) exempted);

typedef struct {
  int32_t target_id;
  int32_t target_length;
  kvec_t(pair_info) reads;
} hiC_target;


void init_hiC_target(hiC_target *ht, int32_t target_id, int32_t target_length);

void free_hiC_target(hiC_target *ht);

typedef struct {
  uint64_t *bits;
  size_t m; // the number of uint64_t holding the bits
  size_t b_n; // the number of bits held
} bit_mask;


// holds two size_ts.. 
typedef struct {
  size_t i;
  size_t n;
}size_t_pair;

/*! @typedef
  @abstract Holds data for all selected read pairs. There should be
            one hiC_target for each target defined in the bam header
	    and the order of these should be identical at the beginning.
	    
	    However, it might be desirable to allow re-ordering and
	    of targets as well as re-allocation of read-pairs if
	    the scaffolding is to be changed.
 */
typedef struct {
  kvec_t(int32_t) target_lengths;
  kvec_t(hiC_target) targets;
  // the target_bits and target_ns are used to allow the rollback of
  // kv_push(pair_info, targets[i];
  bit_mask target_bits;
  kvec_t(size_t_pair) target_n;
  // note that alignment pairs must be stored one after the other.
  // with read one before read 2. Which is left and which is right is defined
  // dynamically by vectors of integers.
  kvec_t(sam_record) alignments;
  kvec_t(char*) query_ids;
} hiC_assembly;

hiC_assembly init_hiC_assembly();

void free_hiC_assembly(hiC_assembly *ha);

/*! @function
  @abstract Obtain all read pairs from a bam file sorted by query name (i.e.
            unsorted). The function will fail if the bam is not sorted by name.

  @param bam_file  Name of a bam file
  @param min_sep   Minimum separation between members of a pair (the internal distance)
  @param min_qual  Minimum mapping quality
  @param min_AS    Minimum alignment score.
  @param max_NM    Maximum number of mismatches in alignment.

  @discussion Notes:
      AS and NM may not always be defined. If not defined, then this will not work as
      currently implemented.
 */
hiC_assembly extract_read_pairs(const char* bam_file, int min_sep, int min_qual,
				int min_AS, size_t max_n, int *error);

hiC_assembly extract_read_pairs_2(const char* bam_file, int min_sep, int min_AS, size_t max_pairs,
				  int merge_max_n, int merge_max_size, int max_al_n,
				  kvi exemptions, int *error);

#endif
