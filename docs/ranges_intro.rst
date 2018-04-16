An introduction to Ranges
=========================

Ranges are intervals with associated data that support range operations (like overlaps and intersections) and other methods that are extremely useful for genomic analyses. The ranges can have an arbitrary number of values associated with them. The pyranges objects are called GRanges as an homage to the original Bioconductor objects.

A GRanges-object (genomic range) have ranges which contain information about what Chromosome (more generally called sequence) and optionally strand a range belongs to. This means that when two GRanges-objects are used for range-operations, the operations are performed on each chromosome separately. Strands can be ignored or used to intersect only on the same strand or opposite strand.

Ranges are built from pandas DataFrames and the information in a Range is stored in a dataframe [#]_. This means the vast Python ecosystem for high-performange scientific computing is available to manipulate the data in GRanges-objects.

.. [#] But under the hood a range object has associated high-performance datastructures to support overlapping operations and slicing.

Examples
~~~~~~~~

.. code-block:: python

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

   # +--------------+---------+-------+---------+----------+
   # | Chromosome   |   Start |   End |   Score | Strand   |
   # |--------------+---------+-------+---------+----------|
   # | chr1         |       4 |     7 |   23.8  | +        |
   # | chr1         |       6 |    11 |    0.13 | -        |
   # | chr2         |       0 |    14 |   42.42 | +        |
   # +--------------+---------+-------+---------+----------+
   # GRanges object with 3 sequences from 2 chromosomes.

   gr1["chr1", 0:5]

   # +--------------+---------+-------+---------+----------+
   # | Chromosome   |   Start |   End |   Score | Strand   |
   # |--------------+---------+-------+---------+----------|
   # | chr1         |       4 |     7 |    23.8 | +        |
   # +--------------+---------+-------+---------+----------+
   # GRanges object with 1 sequences from 1 chromosomes.

   gr1["chr1", "-", 6:100]

   # +--------------+---------+-------+---------+----------+
   # | Chromosome   |   Start |   End |   Score | Strand   |
   # |--------------+---------+-------+---------+----------|
   # | chr1         |       6 |    11 |    0.13 | -        |
   # +--------------+---------+-------+---------+----------+
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

   gr2

   # +--------------+---------+-------+---------+----------+
   # | Chromosome   |   Start |   End |   Score | Strand   |
   # |--------------+---------+-------+---------+----------|
   # | chr1         |       5 |     6 |   -0.01 | -        |
   # | chr1         |       9 |    12 |  200    | +        |
   # | chr3         |       0 |    14 |   21.21 | -        |
   # +--------------+---------+-------+---------+----------+
   # GRanges object with 3 sequences from 2 chromosomes.

   gr1.intersection(gr2, strandedness="opposite")

   # +--------------+---------+-------+---------+----------+
   # | Chromosome   |   Start |   End |   Score | Strand   |
   # |--------------+---------+-------+---------+----------|
   # | chr1         |       4 |     7 |   23.8  | +        |
   # | chr1         |       6 |    11 |    0.13 | -        |
   # +--------------+---------+-------+---------+----------+
   # GRanges object with 2 sequences from 1 chromosomes.

   gr1.intersection(gr2, strandedness=False, invert=True)

   # +--------------+---------+-------+---------+----------+
   # | Chromosome   |   Start |   End |   Score | Strand   |
   # |--------------+---------+-------+---------+----------|
   # | chr2         |       0 |    14 |   42.42 | +        |
   # +--------------+---------+-------+---------+----------+
   # GRanges object with 1 sequences from 1 chromosomes.


The range objects also contain other convenience functions.

.. code-block:: python

   gr1.tile(tile_size=5)

   # +--------------+---------+-------+---------+----------+
   # | Chromosome   | Start   | End   | Score   | Strand   |
   # |--------------+---------+-------+---------+----------|
   # | chr1         | 0       | 4     | 23.8    | +        |
   # | chr1         | 5       | 9     | 23.8    | +        |
   # | chr1         | 5       | 9     | 0.13    | -        |
   # | ...          | ...     | ...   | ...     | ...      |
   # | chr2         | 0       | 4     | 42.42   | +        |
   # | chr2         | 5       | 9     | 42.42   | +        |
   # | chr2         | 10      | 14    | 42.42   | +        |
   # +--------------+---------+-------+---------+----------+
   # GRanges object with 7 sequences from 2 chromosomes.

   gr1.cluster()

   # +--------------+---------+-------+---------+----------+-------------+
   # | Chromosome   |   Start |   End |   Score | Strand   |   ClusterID |
   # |--------------+---------+-------+---------+----------+-------------|
   # | chr1         |       4 |     7 |   23.8  | +        |           1 |
   # | chr1         |       6 |    11 |    0.13 | -        |           1 |
   # | chr2         |       0 |    14 |   42.42 | +        |           2 |
   # +--------------+---------+-------+---------+----------+-------------+
   # GRanges object with 3 sequences from 2 chromosomes.
