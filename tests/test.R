source('../read_bam.R')

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
tl <- target.lengths(bam.ptr)

system.time(
    bam.reg <- aligned.region(region, c(0, tl[region]), bam.ptr,
                              opt.flag=2^5-1, transpose=TRUE, ref.seq=ref.seq)
)
## on my laptop:
##  user  system elapsed 
## 0.025   0.000   0.025

## for:
nrow(bam.reg$diff)
## [1] 6594
nrow(bam.reg$al)
## [1] 6494
nrow(bam.reg$ops)
## [1] 15935

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

bam.sc23.f <- "scaffold_23.bam"
sc23 <- load.bam( bam.sc23.f, paste0(bam.sc23.f, ".bam") )
sc23.tl <- target.lengths( sc23 )

system.time(
    sc23.al <- aligned.region("scaffold_23", c(0, sc23.tl["scaffold_23"]), sc23, transpose=TRUE,
                              opt.flag=(1 + 4 + 8 + 5))
)
## on my laptop:
##  user  system elapsed 
## 0.699   0.057   0.756 


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
## bam <- load.bam("rDNA_r5_sub.bam") 
## this includes forward AND reverse reads.
bam <- load.bam("rDNA_r5.bam")
rflag <- 0x1 + 0x2 + 0x4 + 0x8 + 0x10 + 0x20 + 0x200 + 0x800 + 0x1000 + 0x4000
##rflag <-  bitwShiftL( 1, 15 ) - 1
tmp.1 <- sam.read.n(bam, 4000, ret=rflag, resize=TRUE, transpose=TRUE)

ref.fa <- readLines("contig_236.fa")
ref.fa <- paste(ref.fa[-1], collapse="")

## with both forward and reverse reads we cannot simply get the reference sequence
## in the same way:
fwd.mm <- as.data.frame(tmp.1$mm[ bitwAnd( tmp.1$flag[ tmp.1$mm[,'al.i']  ], 16 ) == 0, ])
rev.mm <- as.data.frame(tmp.1$mm[ bitwAnd( tmp.1$flag[ tmp.1$mm[,'al.i']  ], 16 ) == 16, ])

## remove -1 values
fwd.mm <- as.data.frame(fwd.mm[ fwd.mm[,'r.pos'] > 0, ])
rev.mm <- as.data.frame(rev.mm[ rev.mm[,'r.pos'] > 0, ])

dim(fwd.mm) ## [1] 1292114       6
dim(rev.mm) ## [1] 1752534       6

fwd.dn <- with(fwd.mm, substring(ref.fa, r.pos, r.pos + 1))
rev.dn <- with(rev.mm, substring(ref.fa, r.pos-1, r.pos))

sum(fwd.dn == "CG") / length(fwd.dn) ## 0.8873954
sum(rev.dn == "CG") / length(rev.dn) ## 0.8567354

## this is actually pretty good ;-)

## The non-CGs are most likely related to sequencing errors. I have not included the
## base qualities here, but that could also be done. Note that the ML value is conditional
## on the base being correct.

## It doesn't look like there is a very strong pattern from below.
## However what I did note is that that for the first wrong reference position
## the query position coincided exactly with a deletion even (D). That is that
## the beginning of the query was the same as the query looked for.
## I will need to do some serious tracing of the data to work out what the problem is.

## this would seem to be convenient:
mm <- data.frame( tmp.1$mm, ls=ifelse( bitwAnd(16, tmp.1$flag[ tmp.1$mm[,'al.i'] ]) == 0, 0, -1))
mm <- mm[ mm$r.pos > 0, ]


## looks like the problems may be related to specific reads. We can ask if that is the case by:
cg.by.al <- as.data.frame(t(sapply( tapply( 1:nrow(mm), mm[,'al.i'], function(i){
    p <- mm$r.pos[i] + mm$ls[i]
    c( m=sum( substring( ref.fa, p, p+1 ) == "CG"), n=length(p) )
}), eval)))

## this does look like there are a number of sequences that are especially problematic.
with(cg.by.al, hist(m/n, breaks=100))

## lets have quick look at the number of estimates by position;
## do separately for the 
n.by.pos <- tapply( mm[,'r.pos'], mm[,'mod'], function(p){
    tmp <- table(p)
    data.frame( pos=as.numeric( names(tmp)), count=as.numeric(tmp) )
})

## And this looks like it should do:
sapply(n.by.pos, dim)
##        104   109
## [1,] 15282 15282
## [2,]     2     2

## This is probably sufficient to define the real CGs from the wrong ones.
with(n.by.pos[[2]], plot(pos, count))
## note that the coverage in the first part of the data is enormous; more than 1500

## note that the pattern arises as a result of the sorting.
hist(log2(n.by.pos[[2]]$count))

## lets look at the likelihoods by position:
l.by.pos <- tapply( 1:nrow(mm), mm[,'mod'], function(i){
    probs <- seq(0, 1, 0.1)
    tmp <- tapply( mm[i,'mod.l'], mm[i,'r.pos'], function(l){
        c(n=length(l), ml=mean(l), quantile(l, probs=probs)) })
    df <- data.frame( pos=as.numeric(names(tmp)), do.call(rbind, tmp))
    colnames(df) <- c("pos", "n", "ml", paste0("p.", 100*probs))
    df
})

par(mfrow=c(2,1))
### THERE IS SOMETHING INTERESTING WITH THE HYDROXY-METHYLATION!
with(l.by.pos[[1]], plot(pos, ml, cex=log10(n/5), pch=19, col=rgb(0,0,0,0.2)))
with(l.by.pos[[2]], plot(pos, ml, cex=log10(n/5), pch=19, col=rgb(0,0,0,0.2)))

par(mfrow=c(5,1))
par(mar=c(2.1, 4.1, 2.1, 2.1))
with(l.by.pos[[2]], plot(pos, ml, cex=log10(n), pch=19, col=rgb(0,0,0,0.2), main="mean L"))
with(l.by.pos[[2]], plot(pos, p.20, cex=log10(n), pch=19, col=rgb(0,0,0,0.2), main="L q20"))
with(l.by.pos[[2]], plot(pos, p.50, cex=log10(n), pch=19, col=rgb(0,0,0,0.2), main="L q50"))
with(l.by.pos[[2]], plot(pos, p.80, cex=log10(n), pch=19, col=rgb(0,0,0,0.2), main="L q80"))
with(l.by.pos[[2]], plot(pos, p.90, cex=log10(n), pch=19, col=rgb(0,0,0,0.2), main="L q90"))

par(mfrow=c(5,1))
par(mar=c(2.1, 4.1, 2.1, 2.1))
b <- l.by.pos[[2]]$n > 100
with(l.by.pos[[2]][b,], plot(pos, ml, cex=log10(n), pch=19, col=rgb(0,0,0,0.2), main="mean L"))
with(l.by.pos[[2]][b,], plot(pos, p.20, cex=log10(n), pch=19, col=rgb(0,0,0,0.2), main="L q20"))
with(l.by.pos[[2]][b,], plot(pos, p.50, cex=log10(n), pch=19, col=rgb(0,0,0,0.2), main="L q50"))
with(l.by.pos[[2]][b,], plot(pos, p.80, cex=log10(n), pch=19, col=rgb(0,0,0,0.2), main="L q80"))
with(l.by.pos[[2]][b,], plot(pos, p.90, cex=log10(n), pch=19, col=rgb(0,0,0,0.2), main="L q90"))

### 
fwd.i <- with(l.by.pos[[2]], which(pos %in% fwd.mm$r.pos & substring(ref.fa, pos, pos + 1) == "CG"))
rev.i <- with(l.by.pos[[2]], which(pos %in% rev.mm$r.pos & substring(ref.fa, pos-1, pos) == "CG"))

fwd.b <- with(l.by.pos[[2]], (pos[fwd.i] + 1) %in% pos[rev.i])
rev.b <- with(l.by.pos[[2]], (pos[rev.i] - 1) %in% pos[fwd.i])

l.by.pos.fr <- cbind( l.by.pos[[2]][fwd.i[fwd.b], ], l.by.pos[[2]][rev.i[rev.b], ] )
colnames(l.by.pos.fr) <- paste0( c(rep("f.", 14), rep("r.", 14)), colnames(l.by.pos.fr))

## select those with at least 50 or 20 reads (note that these are all CG anyway, so we are ok..
b <- with(l.by.pos.fr, f.n >= 50 & r.n >= 20) ## that is 95 % of all the data anyway.
l.by.pos.fr <- l.by.pos.fr[b,]

par(mfrow=c(2,1))
## this selection says only give me sites where at least some of the CGs have more than
## 50% likelihood of being methylated
b <- with(l.by.pos.fr, f.p.100 > 128 | r.p.100 > 128)
with(l.by.pos.fr[b,], plot( c(f.pos, r.pos), c(f.ml, r.ml), col=c(rep(1, length(f.pos)), rep(2,length(f.pos)))))
with(l.by.pos.fr[b,], segments( f.pos, f.ml, r.pos, r.ml ))
with(l.by.pos.fr[b,], plot( f.pos, log2(f.ml / r.ml) ))
abline(h=0, col='red')


## to get the most likely number of methylated positions
## given a set of likelihoods we can do:
## we can make use of the Poisson Binomial distribution
## (see https://en.wikipedia.org/wiki/Poisson_binomial_distribution)
## there is a package that implements this efficiently (as this is not
## trivial:
## The automatic install failed; due to compilation issues 
## Rcpp: requiring a newer compiler; (source /opt/rh/devtoolset-9/enable)
## PoissonBinomial: requiring fftw3, which I needed to install from source
##                  note that ./configure --enable-shared needs to be used
##                  to create the required .so files
## install.packages("PoissonBinomial")
require("PoissonBinomial")

## obtain forward and reverse likelihoods by position:
## for the two different modifications:
fwd.bp.l <- tapply(1:nrow(fwd.mm), fwd.mm$mod, function(i){ with(fwd.mm[i,], tapply(mod.l, r.pos, eval)) })
rev.bp.l <- tapply(1:nrow(rev.mm), rev.mm$mod, function(i){ with(rev.mm[i,], tapply(mod.l, r.pos, eval)) })
## we are actually only really interested in a subset of these:
## we can take those that are defined in l.by.pos.fr

for(i in 1:length(fwd.bp.l)){
    fwd.bp.l[[i]] <- fwd.bp.l[[i]][ as.numeric(names(fwd.bp.l[[i]])) %in% l.by.pos.fr$f.pos ]
    rev.bp.l[[i]] <- rev.bp.l[[i]][ as.numeric(names(rev.bp.l[[i]])) %in% l.by.pos.fr$r.pos ]
}

length(fwd.bp.l[[1]]) ##1243
length(fwd.bp.l[[2]]) ##1243
length(rev.bp.l[[1]]) ##1243

## tapply does do what we want with integers as levels
all(names(fwd.bp.l[[2]]) == as.character(l.by.pos.fr$f.pos)) ## TRUE
all(names(rev.bp.l[[2]]) == as.character(l.by.pos.fr$r.pos)) ## TRUE

which.max(dpbinom( NULL, (0.5 + fwd.bp.l[[1]]) / 256, method="DivideFFT" ))
## this gives 46, suggesting 45 is the most likely. Note that almost all the available
## methods: "DivideFFT", "Convolve", "Characteristic", "Recursive", "Mean", "GeoMean",
##          "GeoMeanCounter", "Poisson", "Normal" or "RefinedNormal"
## also give 44. Notable exceptions are: GeoMean: 14
##                                       GeoMeanCounter: 130
## the shape of the likelihoods is similar for DivideFFT, Convolve and Characteristic
## Mean, and Poisson give less sharp peaks.

## how fast is this then:
## this takes a few seconds for each one. Note that although the number of positions is not large,
## the sequencing depth is unusually so, and this may slow down the process more.
## however, we should also provide a confidence range;

## x is a set of probabilities or likelihoods.
n.met <- function(x, cnf=0.95){
    ## note that I do not add 0.5 to x here; this is because it
    ## is common to have a large number of low likelihoods of methylation
    ## that I have a suspicion will inflate the estimate of total numbers
    ## of methylated positions.
    l <- dpbinom(NULL, x/256, method="DivideFFT")
    o <- order(l, decreasing=TRUE)
    cnf.i <- which( cumsum(l[o]) >= cnf )[1]
    tmp <- c( ml=o[1], l=min(o[1:cnf.i]), h=max(o[1:cnf.i]) )
    c(n=length(x), tmp - 1)
}

fwd.bp.ml <- lapply(fwd.bp.l, function(x){ as.data.frame(t(sapply(x, n.met))) })
rev.bp.ml <- lapply(rev.bp.l, function(x){ as.data.frame(t(sapply(x, n.met))) })

## To see the relationship between uncertainty, position and methylation
## level.
with(fwd.bp.ml[[1]], plot(l.by.pos.fr$f.pos, ((h-l)/n) / (ml/n), cex=4 * ml/n ))
with(rev.bp.ml[[1]], points(l.by.pos.fr$r.pos, ((h-l)/n) / (ml/n), cex=4 * ml/n, col='red' ))
##
with(fwd.bp.ml[[2]], plot(l.by.pos.fr$f.pos, ((h-l)/n) / (ml/n), cex=4 * ml/n ))
with(rev.bp.ml[[2]], points(l.by.pos.fr$r.pos, ((h-l)/n) / (ml/n), cex=4 * ml/n, col='red' ))

## fwd.bp.ml <- as.data.frame(t(sapply(fwd.bp.l, function(x){ c(n=length(x), m=which.max(dpbinom(NULL, (0.5 + x)/256, method="DivideFFT"))-1) })))
## rev.bp.ml <- as.data.frame(t(sapply(rev.bp.l, function(x){ c(n=length(x), m=which.max(dpbinom(NULL, (0.5 + x)/256, method="DivideFFT"))-1) })))

par(mfrow=c(2,1))
with(l.by.pos.fr, plot( f.pos, log2(f.ml / r.ml), main="log2 Mean likelihood ratio" ))
abline(h=0, col='red')
## 
plot(l.by.pos.fr$f.pos, log2((fwd.bp.ml[[2]]$m / fwd.bp.ml[[2]]$n) / (rev.bp.ml[[2]]$m / rev.bp.ml[[2]]$n)),
     main="log2 ML methylation rate ratio",
     cex=2 * pmax(fwd.bp.ml[[2]]$m / fwd.bp.ml[[2]]$n, rev.bp.ml[[2]]$m / rev.bp.ml[[2]]$n))
abline(h=0, col='red')
legend('topleft', legend=c(0.25, 0.5, 0.75, 1), pch=1, pt.cex=2 * c(0.25, 0.5, 0.75, 1))


## we can also try the following to see the numbers are distributed at individual sites
k <- 1
o <- with(fwd.bp.ml[[k]], order(ml/n, decreasing=TRUE))
o <- with(rev.bp.ml[[k]], order(ml/n, decreasing=TRUE))
for(i in o){
    par(mfrow=c(2,1))
    par(mar=c(5.1, 4.1, 4.1, 4.1))
    fl <- fwd.bp.l[[k]][[i]]
    rl <- rev.bp.l[[k]][[i]]
    fl.l <- dpbinom(NULL, (0.5 + fl)/256, method="DivideFFT")
    rl.l <- dpbinom(NULL, (0.5 + rl)/256, method="DivideFFT")
    ## estimate the 95% confidence range:
    fl.l.o <- order(fl.l, decreasing=TRUE)
    rl.l.o <- order(rl.l, decreasing=TRUE)
    f.95 <- which(cumsum(fl.l[ fl.l.o ]) >= 0.95)[1]
    r.95 <- which(cumsum(rl.l[ rl.l.o ]) >= 0.95)[1]
    plot( sort(fwd.bp.l[[k]][[i]], decreasing=TRUE), main=paste("forward", fwd.bp.ml[[k]]$m[i], rownames(fwd.bp.ml[[k]])[i] ))
    abline(v=fwd.bp.ml[[k]][i,'ml'])
    abline(h=128)
    with(par(), plot.window(xlim=usr[1:2], ylim=range(fl.l), xaxs='i'))
    points( 1:length(fl.l) - 1, fl.l, col='red' )
    abline(v=range( fl.l.o[ 1:f.95 ]), col='blue')
    axis(4)
    plot( sort(rev.bp.l[[k]][[i]], decreasing=TRUE), main=paste("reverse", rev.bp.ml[[k]]$m[i], rownames(rev.bp.ml[[k]])[i] ))
    abline(v=rev.bp.ml[[k]][i,'m'])
    abline(h=128)
    with(par(), plot.window(xlim=usr[1:2], ylim=range(rl.l), xaxs='i'))
    points( 1:length(rl.l) - 1, rl.l, col='red' )
    abline(v=range( rl.l.o[ 1:r.95 ]), col='blue')
    axis(4)
    inpt <- readline(paste("i: ", i))
}


par(mfcol=c(2,3))
plot( with( fwd.bp.ml[[1]], ml/n ), with( rev.bp.ml[[1]], ml/n ), xlab="fwd ho-met", ylab="rev ho-met" )
plot( with( fwd.bp.ml[[2]], ml/n ), with( rev.bp.ml[[2]], ml/n ), xlab="fwd met", ylab="rev met" )
plot( with( fwd.bp.ml[[2]], ml/n ), with( fwd.bp.ml[[1]], ml/n ), xlab="fwd met", ylab="fwd ho-met" )
plot( with( rev.bp.ml[[2]], ml/n ), with( rev.bp.ml[[1]], ml/n ), xlab="rev met", ylab="rev ho-met" )
plot( with( fwd.bp.ml[[2]], ml/n ), with( rev.bp.ml[[1]], ml/n ), xlab="fwd met", ylab="rev ho-met" )
plot( with( rev.bp.ml[[2]], ml/n ), with( fwd.bp.ml[[1]], ml/n ), xlab="rev met", ylab="fwd ho-met" )


## visualise specified modifications for a specific region and alignments subset.
## mm should be a data.frame like mm that contains al.i, q.pos, mod, mod.n, mod.l, r.pos, ls
## ls indicates the left shift required from the position to give the CG dinucleotide in
## the reference sequence; if -1, it means the read is reversed and that the actual position
## is offset by one.
vis.mod <- function(bam, mm, als, pos, mods, beg, end, border=NA, ref.seq=NULL){
    ## subset als, to only include ones that span the region
    als.span <- cbind('i'=als, 'beg'=bam$pos[als], 'end'=bam$ops[ bam$q.inf[als, 'ops.end'], 'r1' ])
    als.b <- als.span[,'beg'] < end & als.span[,'end'] > beg
    als.span <- als.span[als.b, ]
    mm <- mm[ mm$al.i %in% als.span[,'i'] & (mm$r.pos + mm$ls) %in% pos & mm$mod %in% mods, ]
    mm.fwd <- (mm$ls == 0)
    mm.i <- as.numeric(as.factor(mm$al.i))
    mm.y <- 1:max(mm.i)
    mm.y <- order(mm.fwd[ match(mm.y, mm.i) ])
    mm.y <- match(unique(mm.i), mm.y)
    rects <- as.data.frame(do.call(rbind, tapply(1:nrow(mm), mm$al.i, function(i){
        do.call(rbind, tapply(i, mm$r.pos[i], function(j){
            ## this will define one square, but which can include information from more than one mods
            col <- c(1, 1, 1)
            for(k in 1:length(j))
                col[k] <- 1 - (mm[j[k],'mod.l'] / 256)
            c('x'=mm[j[1], 'r.pos'], 'y'=mm.y[mm.i[i[1]]], 'g'=col[1], 'b'=col[2], 'r'=col[3], 'fwd'=mm.fwd[j[1]])
        }))
    })))
    plot.new()
    plot.window(xlim=c(beg, end), ylim=range(rects[,'y']))
    with(par(), rect(usr[1], usr[3], usr[2], usr[4], col=rgb(0.8, 0.8, 0.8)))
    with(rects, rect( x-0.5, y-0.5, x+0.5, y+0.5, col=rgb(r, g, b), border=rgb(r,g,b)))
    if(!is.null(ref.seq)){
        n.cex=1;
        while(strwidth("A", cex=n.cex) > 0.8)
            n.cex = n.cex * 0.8
        nucs <- substring(ref.seq, beg:end, beg:end)
        with(par(), text(beg:end, usr[3], nucs, pos=3, cex=n.cex))
        with(par(), text(beg:end, usr[4], nucs, pos=1, cex=n.cex))
    }
    axis(1)
    axis(2)
}

## I'm interested in the region between 14,500 to 16,000
## as this contains intersting hydroxymethylation.
ho.beg <- 15000
ho.end <- 16000
##ho.end <- 15000
ho.b <- tmp.1$pos <= ho.beg & tmp.1$ops[ tmp.1$q.inf[,'ops.end'], 'r1' ] >= ho.end
## restrict to very long sequences that map into the end of the sequence:
ho.b <- tmp.1$pos <= 10000 & tmp.1$ops[ tmp.1$q.inf[,'ops.end'], 'r1' ] >= 20000
sum(ho.b) ## 520

al.sel <- sample(which(ho.b), size=200)
pos.sel <- l.by.pos.fr$f.pos
mod.sel <- c(104, 109) ## ho-met, met

par(mfrow=c(1,1))
vis.mod(tmp.1, mm, al.sel, pos.sel, mod.sel[1], ho.beg, ho.end, ref.seq=ref.fa)
vis.mod(tmp.1, mm, al.sel, pos.sel, mod.sel[2], ho.beg, ho.end, ref.seq=ref.fa)

## I think that the putative origin of replication may be here
m.beg <- 10750
m.end <- 11750
## these have to contain the region
m.b.1 <- tmp.1$pos <= m.beg & tmp.1$ops[ tmp.1$q.inf[,'ops.end'], 'r1' ] >= m.end
sum(m.b.1) ## 2722
m.b.2 <- tmp.1$pos < m.end & tmp.1$ops[ tmp.1$q.inf[,'ops.end'], 'r1' ] > m.beg
sum(m.b.2) ## 3737

al.sel <- sample(which(m.b.1), size=200)
par(mfrow=c(1,1))
vis.mod(tmp.1, mm, al.sel, pos.sel, mod.sel[1], m.beg-500, m.end+500, ref.seq=ref.fa)
vis.mod(tmp.1, mm, al.sel, pos.sel, mod.sel[2], m.beg-500, m.end+500, ref.seq=ref.fa)

## lets ask instead..
m.b.1 <- with(tmp.1, pos > (10750-2000) & (bitwAnd(flag, 16) == 0) & ops[ q.inf[,'ops.end'], 'r1'] < (11750+2000))
m.b.2 <- with(tmp.1, (bitwAnd(tmp.1$flag, 16) == 16) & ops[ q.inf[,'ops.end'], 'r1'] < (11750+2000) & pos > (10750-2000))
sum(m.b.1) ## 510
sum(m.b.2) ## 455
sum(m.b.1 | m.b.2) ## 965

al.sel <- sample( which(m.b.1 | m.b.2), size=965 )
vis.mod(tmp.1, mm, al.sel, pos.sel, mod.sel[1:2], m.beg-2000, m.end+2000, ref.seq=ref.fa)
vis.mod(tmp.1, mm, al.sel, pos.sel, mod.sel[2], m.beg-2000, m.end+2000, ref.seq=ref.fa)

plot.met <- function(i){
    cols <- 1:length(i)
    names(cols) <- as.character(i)
    ops <- as.data.frame(with(tmp.1, ops[ ops[,'al.i'] %in% i, ]))
    mm <-  as.data.frame(with(tmp.1, mm[ mm[,'al.i'] %in% i, ]))
    par(mfrow=c(1,1))
    plot.new()
    plot.window(xlim=c(0,255), ylim=range(ops[,c('q0', 'q1')]))
    axis(3)
##
##    with(mm, segments(0, q.pos, mod.l, q.pos, col=cols[as.character(al.i)]))
    with(mm, segments(0, q.pos, mod.l, q.pos, col=rgb(0,0,0,0.1)))
    plot.window(xlim=range(ops[,c('r0', 'r1')]), ylim=range(ops[,c('q0', 'q1')]))
    with(ops, segments(r0, q0, r1, q1, col=cols[as.character(al.i)], lwd=0.5))
    axis(1)
    axis(2)
}
## but really necessary to translate the coordinates.. This is too tricky otherwise.
plot.met( which(m.lp > 0.25) )

plot.met( 400:600 )
