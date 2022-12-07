# read_bam

A few functions to read bam data into R sessions. This probably
overlaps with functionality in rbamtools and rsamtools, but there are
some specific things that I want to play around with starting from a
simple code base. And mainly I want to avoid reading too much
documentation.

To start with the functions will be limited to handling a single bam
indexed bam file, but this can obviously be combined from R.

The code relies on htslib which needs to be installed somewhere on the
target system. Modify `Makevars` to include the path to the htslib
headers, and make sure that the linker is configured to look for
libraries (shared objects, `.so` files) where the htslib coexists (if
root, modify `/etc/ld.so.conf.d` appropriately. If not root, modifying
`LD_LIBRARY_PATH` may also work.  Note though that you should be
familiar with the consequences of changing `LD_LIBRARY_PATH`.

For documentation see the source `read_bam.R`, or failing that
`read_bam.c`. Note that the functions only handle `bam` files indexed
with `bai` indices.

My motivation for writing this is to a large extent to have a *simple*
code base with minimal dependancies that can reasonably efficiently
convert bam data into data structures that allow visualisation and
statistical analyses in `R`. `read_bam.c`, is currently 525 lines of
code, which is a little more than I would like, but it is not *too*
difficult to read. Unfortunately the amount of `C` code will probably
grow.

For [very messy] examples, see the `test.R` file.
