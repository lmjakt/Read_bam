## this performs some magic that should not be necessary if created
## as a package
dyn.load( paste(dirname(sys.frame(1)$ofile), "read_bam.so", sep="/") )

load.bam <- function(bam, bam.index=""){
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

clear.bam.iterator <- function(bam.ptr){
    .Call("clear_iterator", bam.ptr)
}

## opt.flag is the bitwise OR of:
## 1 return seq_data as a character vector array
## 2 return positions that differ from reference; requires that ref.seq is specified
## 3 calculate sequencing depth throughout the region
## 4 construct a cigar string.. 
## flag.filter: a vector of two elements, [ required flags, banned flags ]
## min.mq = minimum mapping quality
## min.ql = minimum query length (as reported in the bam_core field; this is not the
##          same as the actual query length.
aligned.region <- function(region, range, bam.ptr, transpose=FALSE, merge=FALSE, flag.filter=c(-1,-1),
                           opt.flag=0L, ref.seq="", min.mq=0, min.ql=0){
    tmp <- .Call("alignments_region", region, as.integer(range), bam.ptr, as.integer(flag.filter),
                 as.integer(opt.flag), ref.seq, as.integer(min.mq), as.integer(min.ql))
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
sam.read.n <- function(bam.ptr, n, ret.f=( 0x1 + 0x200 + 0x400 + 0x800),
                       sel.flags=c(0, 0, 0), resize=FALSE, transpose=FALSE){
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
    r <- intToUtf8( bitwShiftR( alt, 8 ))
    q <- intToUtf8( bitwAnd( alt, 0x000000FF ) );
    matrix(sapply( strsplit(c(r, q), ""), eval), ncol=2)
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



