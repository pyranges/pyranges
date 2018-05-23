An introduction to Rles
=======================

Rles are runlengths - a datastructure that can compactly describe the coverage of a genome (or some other score that is associated with each nucleotide). It is also used for efficient genomewide arithmetic operations on these scores.

Examples
~~~~~~~~

.. runblock:: pycon

   >>> import pyranges as pr
   >>> from pyranges import Rle

   >>> runs = [10, 10, 10, 10]
   >>> values = [0, 1, 0, 0]

   >>> r1 = Rle(runs, values)
   >>> r1

   >>> runs2 = [11, 9, 20]
   >>> values2 = [100, 0, 100]

   >>> r2 = Rle(runs2, values2)
   >>> r2

   >>> r1 + r2

   >>> r1 * r2

   >>> r1.runs
   >>> r1.values

   >>> r1.values += 5

   >>> r2 / r1
