Ranges
======

Ranges are interval datastructures that support range operations (like overlaps and intersections) and other methods that are extremely useful for genomic analyses. The ranges can have an arbitrary number of values associated with them. The pyranges library contains two kinds of ranges, named IRanges and GRanges as an homage to the original Bioconductor objects.

An IRange (integer range) only contains integer ranges, while a GRange (genomic range) have ranges which also contain information about what Chromosome (more generally called sequence) a range belongs to. This means that when two GRanges-objects are used for range-operations, the operations are performed on each chromosome separately.

Ranges are built from pandas DataFrames and the information in a Range is stored in a dataframe [#]_. This means the vast Python ecosystem for high-performange scientific computing is available to manipulate Range-objects.

.. [#] Under the hood a range object has associated high-performance datastructures to support overlapping operations and slicing.

Examples
~~~~~~~~

import pyranges as pr
from pyranges import GRanges

import pandas as pd

from io import StringIO

f1 = """Chromosome Start End Score Strand
chr1 4 7 23.8 +
chr1 6 11 0.13 -
chr2 0 14 42.42 +"""

df1 = pd.read_table(StringIO(f1), sep="\s+")

gr1 = GRanges(df1)

gr1

# +----+--------------+---------+-------+---------+
# |    | Chromosome   |   Start |   End |   Score |
# |----+--------------+---------+-------+---------|
# |  0 | chr1         |       4 |     7 |   23.8  |
# |  1 | chr1         |       6 |    11 |    0.13 |
# |  2 | chr2         |       0 |    14 |   42.42 |
# +----+--------------+---------+-------+---------+
# GRanges object with 3 sequences from 2 chromosomes.

gr1["chr1", 0:5]

# +----+--------------+---------+-------+---------+
# |    | Chromosome   |   Start |   End |   Score |
# |----+--------------+---------+-------+---------|
# |  0 | chr1         |       4 |     7 |    23.8 |
# +----+--------------+---------+-------+---------+
# GRanges object with 1 sequences from 1 chromosomes.

gr1.Score

# 0    23.80
# 1     0.13
# 2    42.42
# Name: Score, dtype: float64

f2 = """Chromosome Start End Score Strand
chr1 5 6 -0.01 -
chr1 9 12 200 +
chr3 0 14 21.21 -"""

df2 = pd.read_table(StringIO(f2), sep="\s+")

gr2 = GRanges(df2)

gr1.intersection(gr2)
