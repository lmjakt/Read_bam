source('read_bam.R')

read.fasta <- function(fn){
    lines <- readLines(fn)
    id.i <- grep("^>", lines)
    s.beg <- id.i + 1
    s.end <- c(id.i - 1, length(lines))[-1]
    seq <- mapply(function(b, e){
        paste(lines[b:e], collapse="")
        }, s.beg, s.end)
    names(seq) <- sub("^>([^ ]+).+", "\\1", lines[id.i])
    seq
}
ref.seq <- read.fasta("KI270733.fasta")

bam.f <- "test_srt.bam"
bam.i <- "test_srt.bam.bai"
region <- "KI270733.1"
range <- c(130000, 131000)

bam.ptr <- load.bam(bam.f, bam.i)
bam.reg <- aligned.region(region, range, bam.ptr, opt.flag=1L, transpose=TRUE)
tmp1 <- aligned.region(region, range, bam.ptr, opt.flag=3L)

tmp <- aligned.region("1", range, bam.ptr, opt.flag=1L)
tmp <- aligned.region('contig_1', range, bam.ptr, opt.flag=1L)

bam.reg <- aligned.region(region, range, bam.ptr, opt.flag=7L, ref.seq=ref.seq, transpose=TRUE)

bam.reg <- aligned.region(region, range, bam.ptr, opt.flag=7L, ref.seq=ref.seq, transpose=TRUE)

bam.reg <- aligned.region(region, range, bam.ptr, opt.flag=3L, transpose=TRUE)


tlengths <- target.lengths(bam.ptr)
tlengths[region]
range.2 <- c(1, tlengths[region])

system.time(
        bam.reg.2 <- aligned.region(region, range.2, bam.ptr, transpose=TRUE, flag.filter=c(16, 0), opt.flag=1L)
    )
## user  system elapsed 
## 0.02    0.00    0.02 

system.time(
        bam.reg.2 <- aligned.region(region, range.2, bam.ptr, transpose=TRUE, flag.filter=c(16, 0), opt.flag=3L, ref.seq=ref.seq)
)

## user  system elapsed 
## 0.021   0.001   0.021 


bam.reg.3 <- aligned.region(region, range.2, bam.ptr, transpose=TRUE, opt.flag=7L, ref.seq=ref.seq)

ars.ref.f <- paste("/home/lmj/genomes/lophius/ont_data/assembly/flye/v2/ref_to_v2_minimap/", c("flye_v2_ref.bam", "flye_v2_ref.bam.bai"), sep="") 
bam.ptr.2 <- load.bam(ars.ref.f[1], ars.ref.f[2])

tmp <- aligned.region("contig_77", c(1, 1000000), bam.ptr.2, transpose=TRUE, opt.flag=7L)

system.time(
        bam.reg.2 <- aligned.region(region, range.2, bam.ptr, transpose=TRUE, flag.filter=c(0, 0), opt.flag=3L, ref.seq=ref.seq)
)
##  user  system elapsed 
## 0.031   0.001   0.032 

system.time(
        bam.reg.2 <- aligned.region(region, range.2, bam.ptr, transpose=TRUE, flag.filter=c(0, 0), opt.flag=7L, ref.seq=ref.seq)
)
##  user  system elapsed 
## 0.031   0.001   0.032 


system.time(
    for(i in 1:500){
        bam.reg.2 <- aligned.region(region, range.2, bam.ptr, transpose=TRUE, flag.filter=c(16, 0), opt.flag=1L)
    }
)

## user  system elapsed 
## 5.970   0.057   6.026 

plot.new()
with(bam.reg, plot.window(xlim=range(ops[c(4,6),]), ylim=range(ops[c(5,7),])))
with(bam.reg, segments( ops[4,], ops[5,], ops[6,], ops[7,], col=ops[2,] ))
axis(1)
axis(2)


plot.new()
with(bam.reg.2, plot.window(xlim=range(ops[,c(4,6)]), ylim=range(ops[,c(5,7)])))
with(bam.reg.2, segments( ops[,4], ops[,5], ops[,6], ops[,7], col=ops[,2] ))
axis(1)
axis(2)



sam.flags()
sam.flags(83)

## to get the operations codes:
cigar.str <- .Call("bam_cigar_str")
nuc.table <- .Call("nuc_table")

i <- tmp$ops['op',]
substring(cigar.str, i, i)


source('read_bam.R')
### try an unindexed bam file..
bam.f <- "test_srt.bam"
bam <- load.bam(bam.f, "test_srt.bam.bai")

tmp <- sam.read.n(bam, 10000, ret.f=(2^12-1))
tlen <- target.lengths( bam )

tlen[ unique(tmp$ref) ]

bam <- set.bam.iterator( bam, "3" )
tmp <- sam.read.n(bam, 10000, ret.f=(2^12-1))
tmp$n
## [1] 903

bam <- clear.bam.iterator(bam)
tmp <- sam.read.n(bam, 10000, ret.f=(2^12-1))
tlen <- target.lengths( bam )


table(tmp$ref)

## /bam.f <- "~/genomes/lophius/ont_data/EBP_v14/BF2_27062023/BF2_27062023_dorado_duplex/BF2_dorado_duplex.bam"
bam.f <- "BF2_dorado_duplex.bam"
file.exists(bam.f)

bam <- load.bam(bam.f, "")

## flag stats:
system.time( fstat <- sam.flag.stats(bam, 3) )
##    user  system elapsed 
## 191.206   4.670 195.867 

bam <- load.bam(bam.f, "")
system.time( fstat <- sam.flag.stats(bam, 1) )
##    user  system elapsed 
## 185.480   4.784 190.245 

bam <- load.bam(bam.f, "")
system.time( fstat <- sam.flag.stats(bam, 2) )
##    user  system elapsed 
## 184.936   4.568 189.486 

tmp <- list(count=1);
total <- 0
system.time(
    while(tmp[[1]] > 0){
        tmp <- sam.read.n(bam, 1000)
        total <- total + tmp[[1]]
    }
)
##    user  system elapsed 
## 304.975   6.864 311.811 
## not so bad, compared to the hours of loading fastq.gz
## that's 5 minutes.

tmp.dx <- lapply(1:5, function(x){
    tmp <- sam.read.n(bam, 1000)
## check the dx tag to see if the reads have been given as duplexes:
    dx <- as.numeric(sub(".*?dx:i:([0-9]).*", "\\1", tmp[[5]]))
    c(tmp, list('dx'=dx))
})

tmp1 <- sam.read.n(bam, 100)
tmp2 <- sam.read.n(bam, 100)


#### test bcf reading:
source('read_bam.R')
bcf <- load.bcf("Astac_LS420026.bcf", "Astac_LS420026.bcf.csi")
tmp <- bcf.read.n(bcf, 10, c("GT", "AD", "DP", "GQ"))


## test of MM aux parsing;
source('read_bam.R')
bam <- load.bam("rDNA_r5.bam")
rflag <- 0x1 + 0x2 + 0x4 + 0x8 + 0x10 + 0x20 + 0x200 + 0x800 + 0x1000 + 0x4000
##rflag <-  bitwShiftL( 1, 15 ) - 1
tmp.1 <- sam.read.n(bam, 1000, ret=rflag)

## for convenience:
tmp.1$mm <- t(tmp.1$mm)
## check the positions
tmp.pos <- with(tmp.1, tapply(1:nrow(mm), mm[,'al.i'], function(i){ substring( seq[ mm[i,'al.i'] ], mm[i,'q.pos'], mm[i,'q.pos'] + 1 ) }))
table( unlist(tmp.pos) )
##     CG 
## 798828 
## These are all CG. That suggests that we have parsed the file OK.
## but:
range(tmp.1$mm[,'mod.l'])
## 0 255 
hist(tmp.1$mm[,'mod.l'])
## the vast majority of CpGs look unmethylated. 

