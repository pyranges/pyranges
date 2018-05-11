

library(GenomicRanges)
library(data.table)

chip = "/mnt/scratch/endrebak/genomes/chip/UCSD.Aorta.H3K27me3.STL003.bed.gz"
background = "/mnt/scratch/endrebak/genomes/chip/UCSD.Aorta.Input.STL002.bed.gz"

nrows = NULL

c = fread("zcat /mnt/scratch/endrebak/genomes/chip/UCSD.Aorta.H3K27me3.STL003.bed.gz | cut -f 1-3,6", header=FALSE, nrows=nrows)
b = fread("zcat /mnt/scratch/endrebak/genomes/chip/UCSD.Aorta.Input.STL002.bed.gz | cut -f 1-3,6", header=FALSE, nrows=nrows)

print("Starting to create")

start.time <- Sys.time()
c_gr = GRanges(seqnames=c$V1, ranges=IRanges(start=c$V2, end=c$V3), strand=c$V4)
end1.time <- Sys.time()
b_gr = GRanges(seqnames=b$V1, ranges=IRanges(start=b$V2, end=b$V3), strand=b$V4)
end2.time <- Sys.time()

i_gr = intersect(c_gr, b_gr, ignore.strand=TRUE)

end.overlap = Sys.time()

print(end1.time - start.time)
print(end2.time - end1.time)
print(end.overlap - end2.time)
