Run length encoding dicts
=========================

Since you need more than one run length to describe a genome with multiple
chromosomes, pyranges has a datastructure called PyRles for collections of Rles.
It can be created from a PyRanges object by invoking the coverage function.

Rledicts support the arithmetic operations +, -, /, and *.

.. runblock:: pycon

   >>> import pyranges as pr

   >>> gr = pr.load_dataset("chipseq")
   >>> gr_bg = pr.load_dataset("chipseq_background")

   >>> cs = gr.coverage()
   >>> cs

   >>> bg = gr_bg.coverage()
   >>> bg

   >>> cs + bg

When using arithmetic operations with a stranded and an unstranded PyRle, the
stranded PyRle is automatically demoted to an unstranded PyRle.

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("chipseq") # ignore
   >>> gr_bg = pr.load_dataset("chipseq_background") # ignore

   >>> cs = gr.coverage() # ignore

   >>> bg_stranded = gr_bg.coverage(stranded=True)
   >>> bg_stranded

   >>> cs + bg_stranded

Like Rles, PyGRles supports arithmetic operations with numbers.

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("chipseq") # ignore
   >>> cs = gr.coverage() # ignore

   >>> (0.67 + cs) * 5
