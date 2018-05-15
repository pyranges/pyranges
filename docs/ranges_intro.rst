An introduction to PyRanges
===========================

PyRanges are intervals with associated data that support range operations (like
overlap and intersection) and other methods that are useful for genomic
analyses. The ranges can have an arbitrary number of values associated with
them.

A PyRanges object is built from a Pandas dataframe. It should at least contain
the columns Chromosome, Start and End. A column called Strand is optional. Any
other columns in the dataframe are treated as metadata.

The dataframe is kept in the PyRanges object. This means the vast Python
ecosystem for high-performange scientific computing is available to manipulate
the data in PyRanges-objects.

All columns except Chromosome, Start, End and Strand can be changed in any way
you please and more metadata-columns can be added by setting it on the PyRanges
object. If you wish to change the Chromosome, Start, End and Strand columns you
should make a copy of the data from the PyRanges object and use it to
instantiate a new PyRanges.

Examples
~~~~~~~~

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
