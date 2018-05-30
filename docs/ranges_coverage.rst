Turning Ranges into run length encodings
========================================

Ranges can be turned into dicts of run length encodings with the coverage function:

.. runblock:: pycon

   >>> import pyranges as pr
   >>> gr = pr.load_dataset("chipseq")
   >>> gr
   >>> gr.coverage()
   >>> gr.coverage(stranded=True)
