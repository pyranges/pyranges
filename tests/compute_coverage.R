#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

options(warn=-1)

f1 = args[1]
rf = args[2]

df1 = read.table(f1, sep="\t", header=TRUE)

suppressMessages(library(GenomicRanges))

result = coverage(GRanges(seqnames = df1$Chromosome, ranges = IRanges(start = df1$Start, end = df1$End)))

df = data.frame(Runs=runLength(result), Values=runValue(result))

print(df)

write.table(df, rf, sep="\t")
