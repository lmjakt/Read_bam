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
3. `load.bam` is unable to work with the `~` expansion. This should be
   considered as a bug since `R` accepts `~` as a short cut for `$HOME`;
   

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
   **done?** Most (all?) literal number definitions have been moved to common.h
2. `sam.h` defines `BAM_CMATCH`, `BAM_CINS`, etc. I unfortunately redefined a load of `CIG_M`,
   etc. constants in `common.h` (these taking 1 based values). This is messy; I should only
   use the `sam.h` ones. **done**
3. In many cases I defined integer arrays (eg. `int av_column[8]` in `alignments_region()`) that I then
   populate with values from individual variables before passing the array as an argument to `push_column()`.
   This arrays are superfluous as I can instead pass, `(int[]){var1, var2, ...}`. Doing that should make the
   code more compact and more readable, and hence less prone to error.  
   **done**
4. `cigar_to_table()` and `alignments_region()` both parse the cigar data independently. This is because
   `alignments_region()` may need to also update depth and check for mismatches at every position. `sam_read_n()`
   instead creates the table first and then uses it to parse `MM` auxiliary tags (since it needs to do this in
   reverse for reverse complemented query sequences). We could make the functions more consistent though by
   having a helper function that parses a single cigar operation and assigns values to pointers.
   **done** Both `sam_read_n()` and `alignments_region` now use `cigar_to_table()` to parse the cigar data
   and to detect mismatches, base modifications and to calculate depth.
5. `alignments_region()` uses the helper function `bam_seq()` to convert bam sequence data to `R` char* objects.
   `bam_seq()` calls `malloc` each time which then has to be freed; it would be better for it to take a pointer
   to a buffer that it can realloc if needed.
   **complex** This seems inefficient, but if the query sequence is returned to the user, then this pointer is
   anyway stored and it's contents copied to the `R` return `SEXP` object. However, if the query sequence is used
   to detect differences, then we do end up with lots of `malloc` / `free` cycles. So it would make sense to
   simply take the address and copy from that to the return data structure *only* if the query sequence will
   be returned. This needs to be handled with a little bit of care though.
6. In `alignments_region()` I call `strlen()` to determine the length of the reference sequence. I think that
   I should be able to simply use `LENGTH()` on the `CHARSXP` object instead and that this should be more
   efficient. **done**
7. `alignments_region` defines a set of FLAG values that are used to determine what to return. These
   are given as the FLAG values (i.e. 1, 2, 4, ...). It would be better to use bitshift values as
   for `sam_read_n`; these can be chosen so that they correspond to the field in the return list;
   checking them requires expression like: `opt_flag & (1 << AR_Q_DIFF)`, but it means I can do:
   `if(opt_flag ^ (1 << AR_Q_DEPTH)) SET_VECTOR_ELT(ret_data, AR_Q_DEPTH, allocVECTOR(...`. Which
   would remove some dangerous values.
   **partly done** The flag values are now defined by their bit position as for `sam_read_n()`; but 
   I'm not yet using them to define the return data structure.
8. `cigar_to_table()` takes too many arguments. It would be better for all of the arguments to be passed
   as part of a `cigar_parse_options` struct (`cig_opt`).

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
   **partly addressed**
2. One or more R files that run the functions.  
   **partly addressed**
3. Saved results from previous test runs that have been checked for correctness. These should
   be used to compare results after code updates.
4. A file containing a summary of the tests.
   
There is an obvious bootstrap problem here, in that the standardised results have to be created
as part of step 2 whenever a new function is initially tested.


