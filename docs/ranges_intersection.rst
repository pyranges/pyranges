Intersecting Ranges
===================

PyRanges objects can be intersected with other PyRanges to find the subset of
the genome that is contained in both. The regular intersection-method finds the
intersection of all combinations of ranges: [#]_

.. runblock:: pycon

   >>> import pyranges as pr
   >>> gr = pr.load_dataset("aorta")
   >>> gr2 = pr.load_dataset("aorta2")
   >>> gr.intersection(gr2)

The set_intersection method clusters the intervals (i.e. merges them into one) before finding the intersection: [#]_

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("aorta") # ignore
   >>> gr2 = pr.load_dataset("aorta2") # ignore
   >>> gr.set_intersection(gr2)

Both methods also take a strandedness option, which can either be :code:`"same"`, :code:`"opposite"` or :code:`False`/:code:`None`

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gr = pr.load_dataset("aorta") # ignore
   >>> gr2 = pr.load_dataset("aorta2") # ignore
   >>> gr.set_intersection(gr2, strandedness="opposite")

.. [#] This is the same behavior as bedtools intersect.
.. [#] This is the same behavior as Bioconductor GenomicRanges intersect.

The intersection method also takes a how argument, which currently accepts the
option :code:`"containment"`, which requires that the intervals in self be
completely within the intervals in other.

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> f1 = pr.load_dataset("f1")
   >>> f1
   >>> f2 = pr.load_dataset("f2")
   >>> f2
   >>> f2.intersection(f1, how="containment")
