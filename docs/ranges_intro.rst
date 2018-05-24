An introduction to PyRanges
===========================

PyRanges are intervals with associated data that support range operations (like
overlap and intersection) and other methods that are useful for genomic
analyses. The ranges can have an arbitrary number of values associated with
them.

The data in a PyRanges object is stored in a pandas dataframe. This means the
vast Python ecosystem for high-performange scientific computing is available to
manipulate the data in PyRanges-objects.

.. runblock:: pycon

   >>> import pyranges as pr
   >>> from pyranges import PyRanges

   >>> import pandas as pd

   >>> from io import StringIO

   >>> f1 = """Chromosome Start End Score Strand
   ... chr1 4 7 23.8 +
   ... chr1 6 11 0.13 -
   ... chr2 0 14 42.42 +"""

   >>> df1 = pd.read_table(StringIO(f1), sep="\s+")

   >>> gr1 = PyRanges(df1)

   >>> gr1
   >>> gr1["chr1", 0:5]

   >>> gr1["chr1", "-", 6:100]

   >>> gr1.Score

   >>> f2 = """Chromosome Start End Score Strand
   >>> chr1 5 6 -0.01 -
   >>> chr1 9 12 200 +
   >>> chr3 0 14 21.21 -"""

   >>> df2 = pd.read_table(StringIO(f2), sep="\s+")

   >>> gr2 = PyRanges(df2)

   >>> gr2

   >>> gr1.intersection(gr2, strandedness="opposite")

   >>> gr1.intersection(gr2, strandedness=False)

   >>> # The range objects also contain other convenience functions.

   >>> gr1.tile(tile_size=5)

   >>> gr1.cluster()

   >>> gr1.df
