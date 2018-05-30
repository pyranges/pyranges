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

The nearest method also takes a variable how, which can take the values
:code:`"nonoverlapping"`, :code:`"next"`, :code:`"previous"`,
:code:`"next_nonoverlapping"`, :code:`"previous_nonoverlapping"`:. The default
is to find the nearest overlapping interval, no matter which direction.

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> f1 = pr.load_dataset("f1")
   >>> f1
   >>> f2 = pr.load_dataset("f2")
   >>> f2
   >>> f2.nearest(f1, strandedness="opposite", how="next")
   >>> f2.nearest(f1, strandedness="opposite", how="next_nonoverlapping")
