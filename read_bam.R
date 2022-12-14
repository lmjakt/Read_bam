## this performs some magic that should not be necessary if created
## as a package
dyn.load( paste(dirname(sys.frame(1)$ofile), "read_bam.so", sep="/") )

load.bam <- function(bam, bam.index){
    .Call("load_bam", bam, bam.index);
}


## opt.flag is the bitwise OR of:
## 1 return seq_data as a character vector array
## 2 return positions that differ from reference; requires that ref.seq is specified
## 3 calculate sequencing depth throughout the region
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
