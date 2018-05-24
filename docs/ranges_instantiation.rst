Instantiating PyRanges
======================

A PyRanges object can be built in two ways: from a Pandas dataframe or using the
constructor with the seqnames, starts and ends (and optionally strands),
individually.

If you instantiate a PyRanges object from a dataframe, the dataframe should at
least contain the columns Chromosome, Start and End. A column called Strand is
optional. Any other columns in the dataframe are treated as metadata.

.. runblock:: pycon

   >>> import pandas as pd
   >>> import pyranges as pr

   >>> chipseq = pr.get_example_path("chipseq.bed")

   >>> df = pd.read_table(chipseq, header=None,
   ...                    names="Chromosome Start End Name Score Strand".split())

   >>> df.head(2)
   >>> df.tail(2)

   >>> pr.PyRanges(df)

The other way to instantiate a PyRanges object is to use the constructor with keywords:

.. runblock:: pycon

   >>> import pandas as pd # ignore
   >>> import pyranges as pr # ignore

   >>> chipseq = pr.get_example_path("chipseq.bed") # ignore

   >>> df = pd.read_table(chipseq, header=None, names="Chromosome Start End Name Score Strand".split()) # ignore

   >>> pr.PyRanges(seqnames=df.Chromosome, starts=df.Start, ends=df.End)

It is possible to make PyRanges objects out of basic Python datatypes:

.. runblock:: pycon

   >>> import pyranges as pr # ignore

   >>> pr.PyRanges(seqnames="chr1", strands="+", starts=[0, 1, 2], ends=(3, 4, 5))

   >>> pr.PyRanges(seqnames="chr1 chr2 chr3".split(), strands="+ - +".split(),
   ...             starts=[0, 1, 2], ends=(3, 4, 5))
