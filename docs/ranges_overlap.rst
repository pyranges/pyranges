Overlapping Ranges
===================

PyRanges objects can be overlapped with other PyRanges to report the intervals
in self that overlap with those in other.

.. runblock:: pycon

   >>> import pyranges as pr
   >>> gr = pr.load_dataset("aorta")
   >>> gr2 = pr.load_dataset("aorta2")
   >>> gr.overlap(gr2)

Both methods also take a strandedness option, which can either be :code:`"same"`, :code:`"opposite"` or :code:`False`/:code:`None`

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("aorta") # ignore
   >>> gr2 = pr.load_dataset("aorta2") # ignore
   >>> gr.overlap(gr2, strandedness="opposite")
