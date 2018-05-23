Run length encoding dicts
=========================

Since you need more than one run length to describe a genome with multiple
chromosomes, pyranges has a datastructure called PyRles to describe a collection
of Rles. It can be created from a PyRanges object by invoking the coverage function.

Examples
~~~~~~~~

.. runblock:: pycon

>>> import pyranges as pr

>>> gr = pr.load_dataset("chipseq")
>>> gr_bg = pr.load_dataset("chipseq_background")

>>> cs = gr.coverage()
>>> cs

>>> bg = gr_bg.coverage(stranded=True)
>>> bg

>>> cs + bg
