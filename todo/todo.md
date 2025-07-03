# read_bam.c

1. I have many instances of literal numbers in the code. These should be replaced
   with constants defined in `sam.h`, or elsewhere.
2. `sam.h` defines `BAM_CMATCH`, `BAM_CINS`, etc. I unfortunately redefined a load of `CIG_M`,
   etc. constants in `common.h` (these taking 1 based values). This is messy; I should only
   use the `sam.h` ones.
3. `sam_read_n` and `alignments_region` seem to return cigar op tables that are either 0 or 1 based.
   This should be standardised to be one or the other.
4. `alignments_region` should be extended to parse auxiliary strings and MM tags. But for this we
   must make sure that the ops tables are compatible. To me it makes sense to move them to be 0 based
   internally.
5. There is a potential bug in how I encode PHRED quality strings; if not offset, they can contain
   0's. This makes it problematic to encode them as `R` character vectors as these may simply be 0
   terminated byte arrays. Of course one could encode them as integer vectors, but that wastes
   lots of memory. Alternatively one could put 4 values into single ints and extract by doing
   a `cbind()` of four bitwise operations as these are vectorised. But that might just
   postpone the waste of memory.
6. It would seem reasonable for the MM / ML parsing to also extract the base call quality value
   for bases as the modification likelihood is conditional on the base quality. One could also,
   if the reference sequence is avaiable ask whether the base is consistent with the reference.
   If it is, then quite possibly it doesn't really matter.
   
