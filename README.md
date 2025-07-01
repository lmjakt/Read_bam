---
title: read_bam
geometry: "left=3cm,right=3cm"
---

# Motivation

`read_bam` provides a few functions to read and parse bam data into native R
data structures. This probably overlaps with functionality in `rbamtools` and
`rsamtools`, but it differs in that it has minimal external dependancies and
all data is returned as core `R` data structures that can be used directly.

The code relies on `htslib` which needs to be installed somewhere on the target
system. I suspect this makes it difficult to install on a Windows based system.


For documentation and example use of the functions provided, see
further down in this file, or use the source (`read_bam.R`, and
`src/read_bam.c`). For [very messy] examples, see the `tests/test.R`
file.

My motivation for writing this was to a large extent to have a
*simple* code base with minimal dependancies that can reasonably
efficiently convert bam data into data structures that allow
visualisation and statistical analyses using core `R`
functions. `read_bam.c`, is currently around 2000 lines of code, which
is a little more than I would like, but it is not *too* difficult to
read. Unfortunately the amount of `C` code will probably
grow. 


Simplicity in this case isn't so much the amount of code, but
primarily a matter of reducing the number of levels of abstraction
used. Reducing this means that it is easier to find the code that
implements a particular function; and this in turn (to me anyway)
makes it easier to modify the code to do different things when
needed. I find this approach easier compared to either using or
defining universally useful APIs.

# Installation (compilation)

To use `read_bam` you will need to:

1. Install the `htslib` library. Instructions are available from the `htslib`
   github repository. You will need to be aware of where the `include` and
   `library` files are installed to. If you are root and use the usual,
   `./configure`, `make` and `make install`, then this will be in
   `/usr/local/`. If you do not have root access, have a look at the output of
   `./configure --help`. You can set the `prefix` to a location in your own
   home directory.
2. Make the `htslib` library available to the linker. To do this, read up on
   the use of the environment variable `LD_LIBRARY_PATH` and the `ldconfig`
   program if you don't know what this means. If you have a working `samtools`
   installation then it is likely that this has already been done. But you need
   to make sure that the header files that are included during installation
   refer to the library that the linker will load.
3. Compile the `read_bam` source code in order to create a shared object
   (`.so`) file (on windows these are known as `.dll` for dynamically linked
   libraries). Fortunately, `R` has a function for doing this and as long as
   you are using a linux system it is generally straightforward.

To compile the `read_bam` library, change to the `src/` subdirectory. Have a
look at the `Makevars` file. On my system it contains the following:

```sh
## CPPFLAGS=-I/usr/local/include
## CFLAGS=-lhts
PKG_LIBS=-lhts
```

The two first lines are commented out; these contain options that
should be passed to `gcc` (the compiler) when compiling the source
code (I should remove these at some point in the future since they
should not be used in `Makevars` files). `CPPFLAGS` contain options
related to the `C` pre-processor; here it specifies a directory that
may contain header files; the second one states that the executable
should be linked to the `libhts.so` library.

If you have installed `libhts` to a non-standard directory, then you
should probably include something like:

```sh
PKG_CPPFLAGS='-I<libhts installation directory>'
```

in the `Makevars` file. I haven't actually tried that, but something like
that should work.

Change to the `src` subdirectory and use `R CMD SHLIB` to create the library
file:

```sh
## make sure you have changed directory to the place you have
## downloaded the files to and:
cd src
R CMD SHLIB read_bam.c common.c
```

This should produce three new files: `read_bam.o`, `common.o` and
`read_bam.so`. The last of these three is the library file. This needs
to be loaded in order to use the compiled functions. This is done
using the `dyn.load`, which is called as the first line from the
`read_bam.R` file.  Hence you can simply source the `read_bam.R` file
in order for the functions to be available (see usage).

# Usage

To use the functions. From an `R` session:

```R
source("<path_to_code/read_bam.R>")
```

You may get some error messages about `libhts` not being available; if this
is the case then it means that your linker has not been configured to look
for `.so` files in the location where `libhts` has been installed.

Sourcing `read_bam.R` defines a number of wrapper and accessory functions into
your workspace. This is different from what happens when you use `library` to
import a package; in particular, all functions that are defined in
`read_bam.R` can be seen as normal `R` objects when you call `ls()`; this is
not always ideal, but I currently have no desire to make the functions
available as part of a proper `R` package as doing so precludes changing the
API in future versions.

To access data in `bam` or `bcf` files you first need to create a handle for
the files; this is done using the `load.bam` or `load.bcf` functions. These
return what are known as external pointers; you cannot do very much with them
yourself, but you will need to pass them as arguments to functions accessing
data from `bam` and `bcf` files. Note that the external pointers do not
survive R save and reload cycles. They will hence need to be recreated if you
save and then reload your image.

## A use example

If I have installed `read_bam` in `~/R/Read_bam`, have a sorted and indexed
bam file in my current working directory called `test.bam` with an index
called `test.bam.bai` then the following can be used to obtain all alignments
to the longest reference sequence:

```R
## import the functions using source
source("~/R/Read_bam/read_bam.R")

## Create the external pointer
bm <- load.bam("test.bam", "test.bam.bai")

## Obtain the lengths of the reference sequences as
## a named integer vector:
ref.lengths <- target.lengths(bm)

## sort ref.lengths long to short
ref.lengths <- sort(ref.lengths, decreasing=TRUE)

## obtain all alignments from the longest reference:
als <- aligned.region(names(ref.lengths)[1], c(0, ref.lengths[1]),
                      transpose=TRUE)
## als, is a named list containing a number of different entries.
## many of which will be NULL in this case. For details see the
## detailed explanations of the function parameters further down.
```

# Functions provided

## `load.bam(<bam_file>, <bamindex_file>)`

Opens a `bam` file and associates it with an external pointer.

Arguments:

1. `bam_file`: the name of a bam file.  
    Required.
2. `bamindex_file`: the name of a bam index file.  
   Optional, if not specified, defaults to the value of `bam_file` with a
   `.bai` suffix appended. If the file does not exist, then functions using
   the index will not be available.

Value:

An external pointer that should be passed to functions that read from bam
files.

Note that the function may also work with `sam` files, though without any
functions requiring an index.

## `load.bcf(<bcf_file>, <bcf_index>)`

Opens a `bcf` file and associates it with an external pointer.

Arguments are as for `load_bam`, except that the name of the index file
defaults to an empty string and must be specified if index operations are
to be used.

## `set.bam.iterator(<bam.ptr>, <region>, <range>`)

Sets the read position of a sorted and indexed bam file to the position
specified by `<region>` and `<range>` (optional).

Arguments:

1. `bam.ptr`: an external pointer as returned by `load.bam()`.
2. `region`: the name of a reference sequence.
3. `range`: two integer values specifying the begin and end of the region from
   which alignments should be obtained.

Value:

An external pointer; it should not be necessary to handle this as the function
changes the state of its first argument.

## `clear.bam.iterator(<bam.ptr>)`

Resets the read position of the external pointer given as its first argument.

## `aligned.region()`

Returns alignment data for a specified region. This requires that the bam file
has been sorted and that an index has been created and specified when creating
the bam handle (`load.bam()`).

### Arguments:

1. `region`: the name of a reference sequence (a single element character
   vector).
2. `range`: the beginning and end positions from which alignments should be
   extracted. Should be given as a two element numeric vector.
3. `bam.ptr`: an external pointer as returned by `load.bam()`.
4. `transpose=TRUE`: If true, transpose matrix data structures to have a fixed
   number of columns and variable numbers of rows.
5. `flag.filter=c(-1,-1)`: A numeric vector of two elements: required and banned
   flags. These act on the sam flag field to include and exclude specific
   types of alignments. Users need to understand bitwise flags in order to
   specify these. The default value accepts alignments with any flag set.
6. `opt.flag=0`: A bitwise flag specifying optional data that can be returned
   in addition to the primary alignment coordinates. The bits used are:


	| bit | decimal | hex | Description |
	| --: | -----:  | --: | :---------- |
	| 1 |  1 | 0x1 | Return query sequences as a character vector named `query` |
    | 2 |  2 | 0x2 | Return positions that differ from a reference sequence supplied by the ref.seq argument. |
    | 3 |  4 | 0x4 | Calculate sequencing depth throughout the specified region. |
    | 4 |  8 | 0x8 | Construct cigar strings for alignments. |
    | 5 | 16 | 0x10 | Return individual base qualities. |

	The argument is constructed by a bitwise `OR` of the individual bits.
7. `ref.seq=""`: The sequence of the region specified as the first
   argument. This sequence is required if the `opt.flag` includes 0x2. It
   should be given as a single element character vector.
8. `min.mq=0`: The minimum mapping quality.
9. `min.ql=0`: The minimum query length.

### Value:

`aligned_region()` returns a named list containing the following elements:

1. `ref`: The name of the region specified. This will always be a single
   element character vector.
2. `query`: The names of the query sequences for the returned alignments. A
   character vector with one entry for every alignment.
3. `al`: A numeric matrix. If `transpose` is `TRUE`, then it will have a
   fixed number of columns:

   | Column | Description |
   | --- | -------------- |
   | flag | The sam flag of the alignment |
   | r.beg | The reference start position |
   | r.end | The reference end position |
   | q.beg | The query start position |
   | q.end | The query end position |
   | mqual | The mapping quality |
   | qlen  | The query length. Note that this is affected by hard clipping and may not represent the read length. |
   | qclen | The query length inferred from the cigar data. This should be the same as qlen if the cigar data does not include any hard clip operations.|


4. `ops`: The cigar data for all alignments encoded as a single integer
   matrix. If `transpose` is `TRUE` it will have the following columns:


   | Column | Description |
   | --- | -------------- |
   | al.i | The index of the alignment for the current operation. If `transpose` is `TRUE`, then this will correspond to the row number in the `al` matrix.|
   | op | The cigar operation represented as an integer value with "MIDNSHP=XB" mapping to the range 1-10.|
   | type | 1, 2, or 3, depending on whether the operation consumes the query, reference or both respectively. |
   | r0 | The reference op start position (1 based) |
   | q0 | The query op start position (1 based) |
   | r1 | The reference op end position plus one (such that r1 - r0 gives the length of the operation if the type is 2 or 3). |
   | q1 | The query op end position plus one. |
   | op.l | The op length |


5. `seq`: If `opt.flag` includes 0x1 then a character vector giving the query
   sequences as present in the bam file. Otherwise `NULL`.

6. `diff`: `NULL` If `opt.flag` does not include 0x2. Otherwise, a matrix
   giving the positions and bases in the reference and query of any
   mis-matched bases. The columns or rows of this matrix are:
   
   | Column | Description |
   | --- | -------------- |
   | al.i | The index of the alignment for the current position. If `transpose` is `TRUE`, then this will correspond to the row number in the `al` matrix. |
   | r.pos | The reference position (1 based). |
   | q.pos | The query position (1 based). |
   | nuc | The reference and query bases and the query quality encoded in the lower 24 bits of a 32 bit integer. These can be extracted using the alt.nuc.q()` function. |
	

7. `depth`: NULL if `opt.flag` does not include 0x4. Otherwise an integer
   vector the same length as the region requested. Each entry in the vector
   gives the sequence depth at the corresponding reference position.

8. `cigar`: NULL if `opt.flag` does not include 0x8. Otherwise a character
   vector with one element for each alignment holding cigar strings for the
   alignments.

9. `qual`: NULL if `opt.flag` does not include 0x10. Otherwise a character
   vector with one element for each alignment. Each element contains a
   (usually) non-printable string. This can be converted to integer quality
   values using the `utf8ToInt()` function. (POTENTIAL BUG: 0 quality values
   may end up being interpreted as end of strings. I need to confirm this).

The key elements of this list are the `al` and `ops` matrices. In particular
the `ops` makes it possible to visualise all alignments with a single call to
`segments` after setting up a suitable plotting surface. For example:

```R
## here als, is a list returned by a call to aligned.region()
## with transpose=TRUE

## set up a plotting surface
plot.new()
with(als, plot.window(xlim=range(al[,c("r.beg", "r.end")]),
                      ylim=c(0, max(al[,"q.end"]))))
## draw all alignments:
with(als, segments(ops[,'r0'], ops[,'q0'], ops[,'r1'], ops[,'q1']))

```

!["A simple visualisation of alignments"](tests/al_plot_01.png "A simple visualisation of alignments")

The `ops` and `diff` matrices can also be used to identify potential variants. To
find locations that have an excess of indels one can simply use the `R` `tapply()`
and `table()` functions:

```R
indel.pos <- with(als, tapply(ops[,'r0'], ops[,'op'], table))
```

Similarly to find locations mismatch excesses (possibly indicating heterozygosity
position) one can:

```R
mm.pos <- table(als$diff[,'r.pos'])
```

And again, this can be combined with table to ask what types of mismatches are
more common:

```R
mm.n.pos <- with(als, tapply(diff[,'r.pos'], bitwAnd(0xFFFF, diff[,'nuc']), table))
```

Here `bitwAnd` of the `nuc` column and `0xFFFF` is used to extract the lowest 16 bits 
of the `nuc` column as these hold the reference and query bases (or residues generally).
You can then use core `R` functions to evaluate the evidence for heterozygosity at individual
sites or to identify mutation trends.

## `sam.read.n`

Returns information about a specified number (`n`) of alignments from the
current read position of a bam handle. This differs from `aligned.region()` in
that it can be used with unsorted and un-indexed bam files. However, as a
consequence the function cannot return the depth of a region since alignments
can come from random locations across the genome. `sam.read.n()` also does not
support the identification of base mismatches; this is because to do so, the
user would need to load the entire reference genome into the current `R`
session.

It's interface differs in that the user can select any combination of bam
information fields to return using a bitwise flag as before. `sam.read.n()`
also supports the parsing of base modification data from auxiliary strings.

### Arguments

1. `bam.ptr`: an external pointer as returned by `load.bam()`.
2. `n`: the number of alignments for which information should be returned.
3. `ret.f`: A flag value that determines what fields are extracted from the
   bam file. It is formed by bitwise `OR` of the following bits:

	| bit | decimal | hex | Description |
	| --: | -----:  | --: | :---------- |
	| 1 |  1 | 0x1 | query id |
    | 2 |  2 | 0x2 | sam flag value |
    | 3 |  4 | 0x4 | reference name |
    | 4 |  8 | 0x8 | reference position |
    | 5 | 16 | 0x10 | mapping quality |
	| 6 | 32 | 0x20 | cigar string |
	| 7 | 64 | 0x40 | reference name for next mate |
	| 8 | 128 | 0x80 | position for next mate |
	| 9 | 256 | 0x100 | template length (can be negative) |
	| 10 | 512 | 0x200 | query sequence as given in sam / bam (may be clipped and reverse complemented) |
	| 11 | 1024 | 0x400 | query base qualites |
	| 12 | 2048 | 0x800 | the auxiliary string |
	| 13 | 4096 | 0x1000 | cigar operations as a matrix; this will be set if bit 15 is set. |
	| 14 | 8192 | 0x2000 | *not used* |
	| 15 | 16384 | 0x4000 | parse base modification (MM) data |
	
	`ret.f` defaults to (2^6 - 1 + 0x200 + 0x400 + 0x800); i.e. query id, flag, reference name
	position, mapping quality, query sequence, query qualities and the auxiliary field.
4. `sel.flags`: A vector of three integers (required flags, banned flags and minimium mapping quality).
5. `resize`: a logical value (`TRUE` or `FALSE`) that determines whether the data should be resized
	if fewer than `n` alignments were returned. Defaults to `FALSE`.
6. `transpose`: whether matrix structures returned should be transposed such that they have
   variable numbers of row and fixed numbers of columns (as for `aligned.region()`). Defaults
   to `TRUE`.

### Value:

`sam.read.n` returns a named list with 16 members corresponding to the 15 bits of the flag
values used in `ret.f`. The 14^th^ element returns a matrix containing information about
the query sequence; this is non-null if the cigar string should be parsed to a matrix. The last
element is a single integer giving the number of alignments returned. This can differ from the
lengths of the individual elements in cases where 
All entries in the returned list whose corresponding bits specified in `ret.f` were 0
will be `NULL`. The elements of the list are:

| name | description |
| :--- | :---------- |
|  id  | Query id (character vector) |
| flag |  ref   pos  mapq cigar ref.m pos.m  tlen   seq  qual   aux   ops
