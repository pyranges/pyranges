library(GenomicRanges)
library(data.table)

chip = "/mnt/scratch/endrebak/genomes/chip/UCSD.Aorta.H3K27me3.STL003.bed.gz"
background = "/mnt/scratch/endrebak/genomes/chip/UCSD.Aorta.Input.STL002.bed.gz"


c = fread("zcat /mnt/scratch/endrebak/genomes/chip/UCSD.Aorta.H3K27me3.STL003.bed.gz | head -100", header=FALSE)
b = fread("zcat /mnt/scratch/endrebak/genomes/chip/UCSD.Aorta.Input.STL002.bed.gz | head -100", header=FALSE)

c_gr = PyRanges(seqnames=c$V1, ranges=IRanges(start=c$V2, end=c$V3), strand=c$V6, mcols=cbind(c$V4, c$V5))
b_gr = PyRanges(seqnames=b$V1, ranges=IRanges(start=b$V2, end=b$V3), strand=b$V6, mcols=cbind(b$V4, b$V5))


c = fread("zcat /mnt/scratch/endrebak/genomes/chip/UCSD.Aorta.H3K27me3.STL003.bed.gz | head -10", header=FALSE)
b = fread("zcat /mnt/scratch/endrebak/genomes/chip/UCSD.Aorta.Input.STL002.bed.gz | head -10", header=FALSE)

c_gr = PyRanges(seqnames=c$V1, ranges=IRanges(start=c$V2, end=c$V3), strand=c$V6, mcols=cbind(c$V4, c$V5))
b_gr = PyRanges(seqnames=b$V1, ranges=IRanges(start=b$V2, end=b$V3), strand=b$V6, mcols=cbind(b$V4, b$V5))
