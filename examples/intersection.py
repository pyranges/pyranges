import ray

# ray.init()

print("done init")

import pandas as pd

c = "/mnt/scratch/endrebak/test_modin/chip.bed"
i = "/mnt/scratch/endrebak/test_modin/input.bed"
nrows = None
dtypes = {"Chromosome": "category", "Strand": "category"}
names = "Chromosome Start End Name Score Strand".split()
index_cols = None

from time import time

start = time()
cdf = pd.read_csv(c, sep="\t", header=None, index_col=index_cols, nrows=nrows, names=names, dtype=dtypes)
end = time()
print("read table 1", end - start)


import pyranges as pr


start = time()
cgr = pr.PyRanges(cdf)
end = time()
print("init", end - start)


start = time()
idf = pd.read_csv(i, sep="\t", header=None, index_col=index_cols, nrows=nrows, names=names, dtype=dtypes)
end = time()
print("read input", end - start)

start = time()
igr = pr.PyRanges(idf)
end = time()
print("init2", end - start)

start = time()
result = cgr.intersection(igr)
end = time()
print("intersection", end - start)
print(result)
