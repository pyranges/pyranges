import pandas as pd
import pyranges as pr
f = "/mnt/scratch/endrebak/pyranges_benchmark/data/download/chip_15000000.bed.gz"
f2 = "/mnt/scratch/endrebak/pyranges_benchmark/data/download/input_15000000.bed.gz"
df = pd.read_table(f, header=None, nrows=None)
print("read first")
df2 = pd.read_table(f2, header=None, nrows=None)
print("read second")
df.columns = "Chromosome Start End Name Score Strand".split()
df2.columns = "Chromosome Start End Name Score Strand".split()
gr = pr.PyRanges(df)
gr2 = pr.PyRanges(df2)

import time

import ray
from pyranges.ray import _intersection
# result = _intersection.remote(cdf, cdf, None)

print("intersecting")
start = time.time()
result = gr.intersection(gr2)
end = time.time()
print(end - start)
# result = _intersection.remote()
print("result", result)
# print(ray.get(result))
# print(gr.dfs)
