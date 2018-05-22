Joining Ranges
===================

You can combine all the intervals that overlap in two PyRanges objects with the join method.

.. runblock:: pycon

   >>> import pyranges as pr
   >>> gr = pr.load_dataset("aorta")
   >>> gr2 = pr.load_dataset("aorta2")
   >>> gr.join(gr2)

Both methods also take a strandedness option, which can either be :code:`"same"`, :code:`"opposite"` or :code:`False`/:code:`None`

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("aorta") # ignore
   >>> gr2 = pr.load_dataset("aorta2") # ignore
   >>> gr.join(gr2, strandedness="opposite")

You can also use
