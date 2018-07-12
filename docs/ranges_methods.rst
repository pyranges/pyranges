Methods for manipulating single PyRanges
========================================

There are several methods for manipulating the contents of a PyRanges.

:code:`cluster` is a mathematical set operation which creates a union of all the intervals in the ranges:

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> f1 = pr.load_dataset("f1")
   >>> f1
   >>> f1.cluster()

:code:`tile` turns each interval into one or more subintervals of tile_size.

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> f1 = pr.load_dataset("f1") # ignore
   >>> f1
   >>> f1.tile(tile_size=2)

:code:`tss` finds the starts of the regions (taking direction of transcription into account).

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> f1 = pr.load_dataset("f1") # ignore
   >>> f1.tss()
   >>> f1.tss(slack=5)

:code:`tes` finds the ends of the regions (taking direction of transcription into account).

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> f1 = pr.load_dataset("f1") # ignore
   >>> f1.tes()
   >>> f1.tes(slack=5)
