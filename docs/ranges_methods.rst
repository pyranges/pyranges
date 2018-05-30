Methods for manipulating single PyRanges
========================================

There are currently two methods for manipulating the contents of a PyRanges.

Cluster is a mathematical set operation which creates a union of all the intervals in the ranges:

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> f1 = pr.load_dataset("f1")
   >>> f1
   >>> f1.cluster()

Tile turns each interval into one or more subintervals of tile_size.

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> f1 = pr.load_dataset("f1")
   >>> f1
   >>> f1.tile(tile_size=2)
