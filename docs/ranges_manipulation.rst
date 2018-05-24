Manipulating the data in PyRanges
=================================

The PyRanges try to be a thin, but extremely useful wrapper around genomic data
contained in pandas dataframes. This dataframe is accessible with the df
attribute of the PyRanges object.

.. runblock:: pycon

   >>> import pyranges as pr
   >>> gr = pr.load_dataset("chipseq")
   >>> gr
   >>> gr.df

To access a column of this dataframe, you can ask for the name directly from the
PyRanges object.

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("chipseq") # ignore
   >>> gr.Start

You can directly insert a column by setting the attribute on the PyRanges object:

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("chipseq") # ignore
   >>> gr.stupid_example = "Hi There!"
   >>> gr
   >>> gr.df.drop("stupid_example", axis=1, inplace=True)
   >>> gr

All columns except Chromosome, Start, End and Strand can be changed in any way
you please and more metadata-columns can be added by setting it on the PyRanges
object. If you wish to change the Chromosome, Start, End and Strand columns you
should make a copy of the data from the PyRanges object and use it to
instantiate a new PyRanges object.

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("chipseq") # ignore
   >>> import pandas as pd
   >>> gr.Name = gr.Chromosome + "_" + pd.Series(range(len(gr))).astype(str)
   >>> gr
