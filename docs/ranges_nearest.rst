Finding the closest intervals
=============================

With the nearest-method, you can search for the feature in other that is nearest
the ones in self.

.. runblock:: pycon

   >>> import pyranges as pr
   >>> gr = pr.load_dataset("chipseq")
   >>> gr2 = pr.load_dataset("chipseq_background")
   >>> gr.nearest(gr2, suffix="_Input")

The nearest method takes a strandedness option, which can either be
:code:`"same"`, :code:`"opposite"` or :code:`False`/:code:`None`

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("chipseq") # ignore
   >>> gr2 = pr.load_dataset("chipseq_background") # ignore
   >>> gr.nearest(gr2, suffix="_Input", strandedness="opposite")
