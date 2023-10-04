


Installation
~~~~~~~~~~~~

The preferred way to install pyranges is through the bioconda channel::

    conda install -c bioconda pyranges

You can also try pip::

    pip install pyranges

PyRanges has some dependencies that are optional. They need to be manually installed if you require their functionality:

.. code-block:: none

    # pyfaidx: fetch sequences from fasta
    pip install pyfaidx

    # ray: multicpu
    pip install -U ray

    # pybigwig: write bigwigs
    pip install pybigwig  # or conda install -c bioconda pybigwig

    # bamread: read bam files
    pip install bamread   # or conda install -c bioconda bamread

    # fisher: fast fisher exact
    pip install fisher    # or conda install -c bioconda fisher


Since these are not required for most  pyranges functionality, they are kept separate to prevent the possibility of the pyranges-install failing due to dependencies that fail installation or conflicting dependencies.
