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

The PyRanges documentation is structured in three parts:

#. The tutorial, on the next page, recommended for all new users
#. The how-to pages, further below, where functionalities are grouped by topic
#. The `API reference <https://pyranges.readthedocs.io/>`_, where all methods are explained in detail

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   tutorial
   
   how_to_book

       

.. * :ref:`search`
.. * :ref:`modindex`
