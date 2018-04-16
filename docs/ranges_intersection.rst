Intersecting Ranges
===================

GRanges objects have an intersection method to find overlaps with other GRanges.
It has several options to control how the intersections are performed.

.. code-block:: python

import pyranges as pr

csgr = pr.load_dataset("chipseq")
bggr = pr.load_dataset("chipseq_background")
