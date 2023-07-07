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


