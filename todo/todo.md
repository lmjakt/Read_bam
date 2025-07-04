# `read_bam.c`

## Potential bugs

1. There is a potential bug in how I encode PHRED quality strings; if not offset, they can contain
   0's. This makes it problematic to encode them as `R` character vectors as these may simply be 0
   terminated byte arrays. Of course one could encode them as integer vectors, but that wastes
   lots of memory. Alternatively one could put 4 values into single ints and extract by doing
   a `cbind()` of four bitwise operations as these are vectorised. But that might just
   postpone the waste of memory.
2. In `cigar_to_table()`, the way that the values of `*q_beg` and `*q_end` are set are inconsistent.
   I don't think that these values are ever actually used in the code and since there doesn't appear
   to be a problem I don't think that this is a bug. However, these values are returned to the `R`
   session if `sam_read_n()` is called with the option to parse the cigar data; the question is what
   they should represent; either where the alignment begins, or whether they should indicate the
   position of the first base in the original query string. That is whether they should be sensitive
   to both hard and soft clipping or just hard clipping. I'm leaning towards the latter, but haven't
   made my mind up yet.

## Features

1. `alignments_region` should be extended to parse auxiliary strings and MM tags. But for this we
   must make sure that the ops tables are compatible. To me it makes sense to move them to be 0 based
   internally.
2. It would seem reasonable for the MM / ML parsing to also extract the base call quality value
   for bases as the modification likelihood is conditional on the base quality. One could also,
   if the reference sequence is available ask whether the base is consistent with the reference.
   If it is, then quite possibly it doesn't really matter.
3. It is possible to set an iterator for an unsorted bam file to the beginning of the file using
   a specific value defined in SAM.h. This should be added to allow resetting the iterator.

## Code improvement

1. I have many instances of literal numbers in the code. These should be replaced
   with constants defined in `sam.h`, or elsewhere.  
   **done?**Numeric bam operations have been replaced.
2. `sam.h` defines `BAM_CMATCH`, `BAM_CINS`, etc. I unfortunately redefined a load of `CIG_M`,
   etc. constants in `common.h` (these taking 1 based values). This is messy; I should only
   use the `sam.h` ones.
3. In many cases I defined integer arrays (eg. `int av_column[8]` in `alignments_region()`) that I then
   populate with values from individual variables before passing the array as an argument to `push_column()`.
   This arrays are superfluous as I can instead pass, `(int[]){var1, var2, ...}`. Doing that should make the
   code more compact and more readable, and hence less prone to error.  
   **partly addressed**
4. `cigar_to_table()` and `alignments_region()` both parse the cigar data independently. This is because
   `alignments_region()` may need to also update depth and check for mismatches at every position. `sam_read_n()`
   instead creates the table first and then uses it to parse `MM` auxiliary tags (since it needs to do this in
   reverse for reverse complemented query sequences). We could make the functions more consistent though by
   having a helper function that parses a single cigar operation and assigns values to pointers.

## Irritating inconsistencies

1. `sam_read_n` and `alignments_region` seem to return cigar op tables that are either 0 or 1 based.
   This should be standardised to be one or the other.  
   **done**

## Testing

The `tests` directory should not just contain a complete hodge-podge of code that I have used to
initially test the functions. Instead it should have:

1. Some suitable data files (i.e. `bam`, `bcf`, `sam`, `cram`) that are used to test functions.
   This is partly in place, but there is no description of what these files are, or how they
   can be used.
2. One or more R files that run the functions.
3. Saved results from previous test runs that have been checked for correctness. These should
   be used to compare results after code updates.
4. A file containing a summary of the tests.
   
There is an obvious bootstrap problem here, in that the standardised results have to be created
as part of step 2 whenever a new function is initially tested.


