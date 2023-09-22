.. pyranges documentation master file, created by
   sphinx-quickstart on Fri Jun 23 11:27:13 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introduction
============

Installation
~~~~~~~~~~~~

The preferred way to install pyranges is through the bioconda channel::

    conda install -c bioconda pyranges

You can also try pip::

    pip install pyranges

PyRanges has some dependencies that are optional. They need to be manually installed if you require their functionality::

    pyfaidx: fetch sequences from fasta # pip install pyfaidx
    ray: multicpu   # pip install -U ray
    pybigwig: write bigwigs # pip install pybigwig
                            # or conda install -c bioconda pybigwig
    bamread: read bam files # pip install bamread
                            # or conda install -c bioconda bamread
    fisher: fast fisher exact # pip install fisher
                              # or conda install -c bioconda fisher


Since these are not needed for 99.9% percent of the pyranges functionality, they are kept separate to prevent the possibility of the pyranges-install failing due to dependencies that fail installation or conflicting dependencies.


Citation
~~~~~~~~

http://dx.doi.org/10.1093/bioinformatics/btz615


Documentation outline
~~~~~~~~~~~~~~~~~~~~~


#. 🚀 :doc:`The tutorial <./tutorial>`, on the next page, recommended for all new users
#. 🚀 :doc:`The how-to pages <./how_to_pages>`, further below, where functionalities are grouped by topic
#. 🚀 The `API reference <./autoapi/index.html>`_, where all methods are explained in detail.
#. 🚀 :doc:`The developer guide <./developer_guide>`, to follow in order to contribute to PyRanges.

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:

   tutorial
   
   how_to_pages

   developer_guide

       

.. * :ref:`search`
.. * :ref:`modindex`
