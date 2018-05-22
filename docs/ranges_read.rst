Reading Ranges from common file formats
=======================================

The pyranges library can create PyRanges from three common file formats, namely
gtf, bed and bam. [#]_

bed
~~~

.. runblock:: pycon

   >>> import pyranges as pr
   >>> bed = pr.get_example_path("aorta.bed")
   >>> bed
   >>> pr.read_bed(bed)

bam
~~~


.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> bam = pr.get_example_path("control.bam")
   >>> bam
   >>> pr.read_bam(bam)

gtf
~~~

.. runblock:: pycon

   >>> import pyranges as pr # ignore
   >>> gtf = pr.get_example_path("ensembl.gtf")
   >>> gtf
   >>> pr.read_gtf(gtf)

.. [#] PyRanges uses the pysam library which requires that the bam file must have an index.
