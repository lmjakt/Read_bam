## this performs some magic that should not be necessary if created
## as a package
if(is.loaded("load_bam"))
    dyn.unload( paste(dirname(sys.frame(1)$ofile), "src/read_bam.so", sep="/") )
dyn.load( paste(dirname(sys.frame(1)$ofile), "src/read_bam.so", sep="/") )

load.bam <- function(bam, bam.index=paste0(bam, ".bai")){
    .Call("load_bam", bam, bam.index);
}

load.bcf <- function(bcf, bcf.index=""){
    .Call("load_bcf", bcf, bcf.index)
}

set.bam.iterator <- function(bam.ptr, region, range=NULL){
    if(is.null(range)){
        tlen = target.lengths(bam.ptr)
        if(is.na(tlen[region[1]]))
            stop(paste(region[1], ":no such region"))
        range=c(1, tlen[region[1]])
    }
    .Call("set_iterator", bam.ptr, region, as.integer(range));
}

## This only frees the memory of the iterator and does
## not cause the reading to start from the beginning.
## To start reading from the beginning of a file we
## can set.bam.iterator using a set of special values
## (see sam.h sam_itr_queryi):
## HTS_IDX_NOCOOR, HTS_IDX_START, HTS_IDX_REST, HTS_IDX_NONE
## for the region part. This needs to be implemented in the read_bam.c.
## This is just a reminder about how to implement this.
clear.bam.iterator <- function(bam.ptr){
    .Call("clear_iterator", bam.ptr)
}

## Note opt.flag=1L causes the index to parse the cigar data in
## order to find the reference end and the query start and end
## positions. It turns out that it doesn't make much difference
## to the speed of the index building. So it's on by default.
build.index <- function(bam.ptr, region, opt.flag=1L){
    .Call("build_query_index", bam.ptr, region, as.integer(opt.flag))
}

query.pos <- function(bam.ptr, query.ids){
    tmp <- .Call("query_positions", bam.ptr, query.ids)
    tmp <- data.frame(query=query.ids[tmp[[1]]], query.i=tmp[[1]], t(tmp[[2]]))
    tmp$target.id <- tmp$target.id + 1
    tmp
}

## opt.flag is the bitwise OR of bits:
## 1 return seq_data as a character vector array
## 2 return positions that differ from reference; requires that ref.seq is specified
## 3 calculate sequencing depth throughout the region
## 4 construct a cigar string..
## 5 return qualities
## 6 mate information.
## 7 parse MM auxiliary data (base modifications)
## 8 intron depth vector
## flag.filter: a vector of two elements, [ required flags, banned flags ]
## min.mq = minimum mapping quality
## min.ql = minimum query length (as reported in the bam_core field; this is not the
##          same as the actual query length.
aligned.region <- function(region, range, bam.ptr, transpose=TRUE, flag.filter=c(0,0),
                           opt.flag=0L, ref.seq="", min.mq=0, min.ql=0, max.intron.length=4096L){
    tmp <- .Call("alignments_region", region, as.integer(range), bam.ptr, as.integer(flag.filter),
                 as.integer(opt.flag), ref.seq, as.integer(min.mq), as.integer(min.ql),
                 as.integer(max.intron.length))
    if(transpose){
        tmp <- lapply(tmp, function(x){
            if( length(dim(x)) == 2 )
                return(t(x))
            x
        })
    }
    ## we could consider to merge query and al, into a data.frame
    ## but there is not that much point
    tmp
}

aligned.region.mt <- function(bam.f, bam.f.index, region, range, transpose=TRUE,
                              nthreads=1L, flag.filter=c(0,0),
                              opt.flag=0L, ref.seq="", min.mq=0, min.ql=0, max.intron.length=4096L){
    tmp <- .Call("alignments_region_mt",
          bam.f, bam.f.index, region, as.integer(range), as.integer(flag.filter),
          as.integer(opt.flag), ref.seq, as.integer(min.mq), as.integer(min.ql),
          as.integer(max.intron.length), as.integer(nthreads))
    if(transpose){
        tmp <- lapply(tmp, function(x){
            if( length(dim(x)) == 2 )
                return(t(x))
            x
        })
    }
    ## we could consider to merge query and al, into a data.frame
    ## but there is not that much point
    tmp
}

count.alignments <- function(region, range, bam.ptr, flag.filter=c(-1,-1), min.mq=0){
    .Call("count_region_alignments",
          region, as.integer(range), bam.ptr,
          as.integer(flag.filter),
          as.integer(min.mq))
}

## there are 11 mandatory fields in the bam file; This function was initially
## written to return only id, seq, qual and aux from non-aligned bam files
## However, it is better to use a flag to define the 12 (11 mandatory + aux)
## columns.
## use the following:
## Bit hxd      dec      field
## 1     1        1      id
## 2     2        2      flag
## 3     4        4      rname (reference name)
## 4     8        8      pos
## 5    10       16      mapq
## 6    20       32      cigar
## 7    40       64      rnext (reference for mate)
## 8    80      128      pnext (pos for mate)
## 9   100      256      tlen (template length can be negative)
## 10  200      512      seq
## 11  400     1024      qual
## 12  800     2048      aux
## the default flag is: id, seq, qual, aux
## 13 1000     4096      cigar_operations
## 15 4000     16384     parse MM aux data. 
## if the 13th bit is set, then the cigar string will be parsed and
## a table containing the operations, and a table containing 
## additional query information will be returned (ops and q.inf)
##
## sel.flags are:
## 1. f: only include reads with all flags set
## 2. F: exclude all reads with any flag set
## 3. min_mapq 
sam.read.n <- function(bam.ptr, n, ret.f=(2^6 - 1 + 0x200 + 0x400 + 0x800),
                       sel.flags=c(0, 0, 0), resize=FALSE, transpose=TRUE){
    tmp <- .Call("sam_read_n", bam.ptr, as.integer(n),
                 as.integer(ret.f), as.integer(sel.flags))
    if(resize){
        for(i in 2:length(tmp) - 1){
            if(!is.null(tmp[[i]]) && is.vector(tmp[[i]]))
                tmp[[i]] <- tmp[[i]][ 1:tmp$n ]
        }
    }
    if(transpose){
        for(i in 1:length(tmp)){
            if(is.matrix(tmp[[i]]))
                tmp[[i]] <- t(tmp[[i]])
        }
    }
    tmp
}

bcf.read.n <- function(bcf.ptr, n, fmt.tags=character(), resize=FALSE){
    bcf <- .Call("bcf_read_n", bcf.ptr, as.integer(n), fmt.tags)
    n <- bcf[[ 1 ]]
    i <- which(!names(bcf) %in% fmt.tags)
    j <- which(names(bcf) %in% fmt.tags)
    smp <- bcf[j];
    bcf <- data.frame(bcf[2:max(i)])
    if(resize)
        bcf <- bcf[1:n,]
    list(bcf=bcf, smp=smp, n=n)
}

## ret flag is a bitwise flag that determines whether
## the counts for the indivdiual bits (01) or the
## counts for all possible flags (10) are returned.
## (11 --> both). There does not seem to be any major
## change in the time taken.
sam.flag.stats <- function(bam.ptr, ret.flag){
    tmp <- .Call("sam_flag_stats", bam.ptr, as.integer(ret.flag));
    names(tmp) <- c('bit', 'flag')
    tmp
}

## returns a named integer of target lengths
target.lengths <- function(bam.ptr){
    .Call("target_lengths", bam.ptr)
}

cigar.str <- function(){
    .Call("bam_cigar_str")
}

nuc.table <- function(){
    .Call("nuc_table")
}

## the default will return all of the flag names
sam.flags <- function(flags=0xFFFF){
    decoded <- lapply( flags, function(x){ (1:12)[ as.logical( bitwAnd(2^(1:11), x) )] })
    strings <- .Call("bam_flag", as.integer(flags))
    lapply( 1:length(decoded), function(i){
        data.frame( bit=decoded[[i]], flag=strsplit(strings[i], ",")[[1]], value=2^(decoded[[i]]-1) )
    })
}

## convert an int to two nucs..
alt.nuc <- function(alt){
    r <- intToUtf8( bitwAnd( bitwShiftR( alt, 8 ), 0x000000FF ) )
    q <- intToUtf8( bitwAnd( alt, 0x000000FF ) );
    matrix(sapply( strsplit(c(r, q), ""), eval), ncol=2)
}

## convert an int to two nucs and a quality in a dataframe:
alt.nuc.q <- function(alt){
    r <- intToUtf8( bitwAnd( bitwShiftR( alt, 8 ), 0x000000FF ) )
    q <- intToUtf8( bitwAnd( alt, 0x000000FF ) );
    qual <- bitwAnd( bitwShiftR( alt, 16 ), 0x000000FF )
    rq <- matrix(sapply( strsplit(c(r, q), ""), eval), ncol=2)
    colnames(rq) <- c('r', 'q')
    data.frame( rq, qual=qual )
}

## extract query insertions from the ops table
## this requires that the query sequences have been included
## or it will cause an error
q.ins <- function(reg, is.t=FALSE, margin=0){
    if(!length(reg$seq))
        stop("Alignments data struct does not contain any query sequences")
    if(is.t)
        ops <- t(reg$ops)
    else
        ops <- reg$ops
    margin=abs(margin)
    b <- ops['op',] == 2
    q.i <- ops['al.i',b]
    q.b <- ops['q0',b] - margin 
    q.e <- ops['q1',b] + margin - 1
    seq.ins <- substring( reg$seq[q.i], q.b, q.e )
    data.frame( al.i=q.i, r.pos=ops['r0',b], q.b=q.b, q.e=q.e-1, seq=seq.ins, stringsAsFactors=FALSE )
}

r.ins <- function(reg, ref.seq, is.t=FALSE, margin=0){
    if(is.t)
        ops <- t(reg$ops)
    else
        ops <- reg$ops
    margin <- abs(margin)
    b <- ops['op',] == 3
    q.i <- ops['al.i',b]
    r.b <- ops['r0',b] - margin 
    r.e <- ops['r1',b] + margin - 1
    seq.ins <- substring( ref.seq, r.b, r.e )
    data.frame( al.i=q.i, r.b=r.b, r.e=r.e, q.pos=ops['q0',b], seq=seq.ins, stringsAsFactors=FALSE )
}

## takes auxiliary strings and converts values to more appropriate values.
## aux: a character vector of auxiliary strings
## tags: the tag values that should be extracted
## types: the data types of those tag values (A, C, S, I, f, Z, H) allowed
extract.aux <- function(aux, tags, types){
    tmp <- .Call("extract_aux_tags", aux, tags, types)
    tmp <- do.call(data.frame, tmp)
    colnames(tmp) <- tags
    tmp
}

## translate query positions to reference positions given
## cigar like operations;
## ops should be a matrix of cigar operations as returned by aligned.region
## but note that only:
## op, r0, q0, and q1 will be used. Other values can be set to 0
## qpos should be a matrix with the following columns:
## query position
## first row in the ops table for each query (1-based)
## last row in the ops table for each query (1-based)
## integer (1 if forward, 0 if reverse)
## the query length
query.to.ref.pos <- function(qpos, ops){
    if(!is.integer(qpos))
        qpos <- matrix(as.integer(qpos), nrow=nrow(qpos), ncol=ncol(qpos))
    if(!is.integer(ops))
        ops <- matrix(as.integer(ops), nrow=nrow(ops), ncol=ncol(ops))
    tmp <- .Call("qpos_to_ref_pos", qpos, t(ops))
    if(!is.null(tmp))
        colnames(tmp) <- c("qpos", "rpos", "op_i")
    tmp
}

## extract the 4 individual bytes from each value of an integer value
## and return as a matrix with a column for each byte:
int.to.bytes <- function(i, names=paste0("b", 1:4)){
    tmp <- sapply(3:0, function(j){
        bitwAnd(0xFF, bitwShiftR(i, j * 8))
    })
    colnames(tmp) <- names[1:4]
    tmp
}

## For drawing:
## stolen from simple_range.R: The interface should change at some point.
## x should be a matrix with two columns
## if no columns with name x0 and x1, then the
## first two columns will be taken to represent x0 and x1
## respectively.
arrange.lines <- function(x){
    if(!is.matrix(x) || !is.numeric(x) || ncol(x) < 2)
        stop("x should be a numeric matrix with at least two columns")
    if(!all(c('x0', 'x1') %in% colnames(x)))
        colnames(x)[1:2] <- c('x0', 'x1')
    if(any( x[,'x0'] >= x[,'x1'] ))
        stop("all x0 values must be smaller than x1")
    if(any(is.na(x[,c('x0', 'x1')])))
       stop("NA values not allowed")
    .Call("arrange_lines", as.double(x[,'x0']), as.double(x[,'x1']))
}
