Subsetting Ranges
=================

There are many ways to subset a PyRanges object. Each returns a new PyRanges object and does not change the old one.

.. runblock:: pycon

   >>> import pyranges as pr

   >>> gr = pr.load_dataset("chipseq")


Chromosome only
~~~~~~~~~~~~~~~

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("chipseq") # ignore
   >>> gr["chrX"]

Chromosome and Strand
~~~~~~~~~~~~~~~~~~~~~

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("chipseq") # ignore
   >>> gr["chrX", "-"]

Chromosome and Slice
~~~~~~~~~~~~~~~~~~~~~

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("chipseq") # ignore
   >>> gr["chrX", 150000000:160000000]

Chromosome, Strand and Slice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("chipseq") # ignore
   >>> gr["chrX", "-", 150000000:160000000]

Slice
~~~~~

Only using slices returns all ranges from all chromosomes and strands within those coordinates.

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("chipseq") # ignore
   >>> gr[0:100000]

Strand
~~~~~~

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("chipseq") # ignore
   >>> gr["+"]

Slice and Strand
~~~~~~~~~~~~~~~~

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("chipseq") # ignore
   >>> gr["+", 0:100000]
