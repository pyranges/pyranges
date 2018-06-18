#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

f1 = args[1]
f2 = args[2]
op = args[3]
rf = args[4]

print(args)

print("We are starting in R! We are starting in R! We are starting in R! We are starting in R! We are starting in R! We are starting in R! We are starting in R! ")


library(S4Vectors)


df1 = read.table(f1, sep="\t", header=TRUE)
df2 = read.table(f2, sep="\t", header=TRUE)

sum1 = sum(df1$Runs)
sum2 = sum(df2$Runs)

print(sum1)
print(sum2)
print(sum2 > sum1)

if (sum1 > sum2){
  row = data.frame(sum1 - sum2, 0)
  colnames(row) = c("Runs", "Values")
  df2 = rbind(df2, row)
} else if (sum2 > sum1){
  row = data.frame(sum2 - sum1, 0)
  colnames(row) = c("Runs", "Values")
  df1 = rbind(df1, row)
}

print(df1)
print(df2)

r1 = Rle(df1$Values, df1$Runs)
r2 = Rle(df2$Values, df2$Runs)

print(r1)
print(r2)

print(op)
f = match.fun(op)
## f = *

result = f(r1, r2)

print(result)

df = data.frame(Runs=runLength(result), Values=runValue(result))

write.table(df, rf, sep="\t")

print("We are done in R! We are done in R! We are done in R! We are done in R! We are done in R! We are done in R! We are done in R! ")
