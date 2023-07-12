How-to-book
===========



Introduction to PyRanges
~~~~~~~~~~~~~~~~~~~~~~~~

This is the PyRanges tutorial. For docs, see: `https://pyranges.readthedocs.io/en/latest/`


.. contents:: Contents of How-to pages
   :depth: 3

PyRanges are collections of intervals that support comparison operations (like overlap and intersect) and other methods that are useful for genomic analyses. The ranges can have an arbitrary number of meta-data fields, i.e. columns associated with them.

The data in PyRanges objects are stored in a pandas dataframe. This means the vast Python ecosystem for high-performance scientific computing is available to manipulate the data in PyRanges objects.




  >>> from pyranges import PyRanges
  >>> import pandas as pd
  >>> from io import StringIO
	
  >>> f1 = """Chromosome Start End Score Strand 
  ... chr1 4 7 23.8 +
  ... chr1 6 11 0.13 -
  ... chr2 0 14 42.42 +"""
	
  >>> df1 = pd.read_csv(StringIO(f1), sep="\s+")
  >>> gr1 = PyRanges(df1)


Now we can subset the PyRange in various ways:


  >>> print(gr1)
  +--------------+-----------+-----------+-------------+--------------+
  | Chromosome   |     Start |       End |       Score | Strand       |
  | (category)   |   (int32) |   (int32) |   (float64) | (category)   |
  |--------------+-----------+-----------+-------------+--------------|
  | chr1         |         4 |         7 |       23.8  | +            |
  | chr1         |         6 |        11 |        0.13 | -            |
  | chr2         |         0 |        14 |       42.42 | +            |
  +--------------+-----------+-----------+-------------+--------------+
  Stranded PyRanges object has 3 rows and 5 columns from 2 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
	
  >>> print(gr1["chr1", 0:5])
  +--------------+-----------+-----------+-------------+--------------+
  | Chromosome   |     Start |       End |       Score | Strand       |
  | (category)   |   (int32) |   (int32) |   (float64) | (category)   |
  |--------------+-----------+-----------+-------------+--------------|
  | chr1         |         4 |         7 |        23.8 | +            |
  +--------------+-----------+-----------+-------------+--------------+
  Stranded PyRanges object has 1 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  
  
  >>> print(gr1["chr1", "-", 6:100])
  +--------------+-----------+-----------+-------------+--------------+
  | Chromosome   |     Start |       End |       Score | Strand       |
  | (category)   |   (int32) |   (int32) |   (float64) | (category)   |
  |--------------+-----------+-----------+-------------+--------------|
  | chr1         |         6 |        11 |        0.13 | -            |
  +--------------+-----------+-----------+-------------+--------------+
  Stranded PyRanges object has 1 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
	
  >>> print(gr1.Score)
  0    23.80
  1     0.13
  2    42.42
  Name: Score, dtype: float64
	
	
And we can perform comparison operations with two PyRanges:

  >>> f2 = """Chromosome Start End Score Strand
  ... chr1 5 6 -0.01 -
  ... chr1 9 12 200 +
  ... chr3 0 14 21.21 -"""
	
  >>> df2 = pd.read_csv(StringIO(f2), sep="\s+")
  >>> gr2 = PyRanges(df2)
  >>> print(gr2)
  +--------------+-----------+-----------+-------------+--------------+
  | Chromosome   |     Start |       End |       Score | Strand       |
  | (category)   |   (int32) |   (int32) |   (float64) | (category)   |
  |--------------+-----------+-----------+-------------+--------------|
  | chr1         |         9 |        12 |      200    | +            |
  | chr1         |         5 |         6 |       -0.01 | -            |
  | chr3         |         0 |        14 |       21.21 | -            |
  +--------------+-----------+-----------+-------------+--------------+
  Stranded PyRanges object has 3 rows and 5 columns from 2 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
	
  >>> print(gr1.intersect(gr2, strandedness="opposite"))
  +--------------+-----------+-----------+-------------+--------------+
  | Chromosome   |     Start |       End |       Score | Strand       |
  | (category)   |   (int32) |   (int32) |   (float64) | (category)   |
  |--------------+-----------+-----------+-------------+--------------|
  | chr1         |         5 |         6 |       23.8  | +            |
  | chr1         |         9 |        11 |        0.13 | -            |
  +--------------+-----------+-----------+-------------+--------------+
  Stranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
	
  >>> print(gr1.intersect(gr2, strandedness=False))
  +--------------+-----------+-----------+-------------+--------------+
  | Chromosome   |     Start |       End |       Score | Strand       |
  | (category)   |   (int32) |   (int32) |   (float64) | (category)   |
  |--------------+-----------+-----------+-------------+--------------|
  | chr1         |         5 |         6 |       23.8  | +            |
  | chr1         |         9 |        11 |        0.13 | -            |
  +--------------+-----------+-----------+-------------+--------------+
  Stranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

There are also convenience methods for single PyRanges:

  >>> print(gr1.merge())
  +--------------+-----------+-----------+--------------+
  | Chromosome   |     Start |       End | Strand       |
  | (category)   |   (int32) |   (int32) | (category)   |
  |--------------+-----------+-----------+--------------|
  | chr1         |         4 |         7 | +            |
  | chr1         |         6 |        11 | -            |
  | chr2         |         0 |        14 | +            |
  +--------------+-----------+-----------+--------------+
  Stranded PyRanges object has 3 rows and 4 columns from 2 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

The underlying dataframe can always be accessed:

  >>> print(gr1.df)
  	Chromosome  Start  End  Score Strand
  0       chr1      4    7  23.80      +
  1       chr1      6   11   0.13      -
  2       chr2      0   14  42.42      +



Loading/Creating PyRanges
~~~~~~~~~~~~~~~~~~~~~~~~~


A PyRanges object can be built in four ways:


#. from a Pandas dataframe
#. using the PyRanges constructor with the chromosomes, starts and ends (and optionally strands), individually.
#. using one of the custom reader functions for genomic data (read_bed, read_bam or read_gtf, read_gff3)
#. from a dict (like the ones produced with to_example)


Using a DataFrame
-----------------


If you instantiate a PyRanges object from a dataframe, it should at least contain the columns Chromosome, Start and End. A column called Strand is optional. Any other columns in the dataframe are treated as metadata.


  >>> import pandas as pd
  >>> import pyranges as pr
  >>> chipseq = pr.get_example_path("chipseq.bed")
  >>> df = pd.read_csv(chipseq, header=None, names="Chromosome Start End Name Score Strand".split(), sep="\t")
  >>> print(df.head(2))
  	Chromosome      Start        End Name  Score Strand
  0       chr8   28510032   28510057   U0      0      -
  1       chr7  107153363  107153388   U0      0      -

  >>> print(df.tail(2))
  	Chromosome      Start        End Name  Score Strand
  9998       chr1  194245558  194245583   U0      0      +
  9999       chr8   57916061   57916086   U0      0      +
	 
  >>> print(pr.PyRanges(df))


	
Using constructor keywords
--------------------------


The other way to instantiate a PyRanges object is to use the constructor with keywords:

  >>> gr = pr.PyRanges(chromosomes=df.Chromosome, starts=df.Start, ends=df.End)
  >>> print(gr)
  +--------------+-----------+-----------+
  | Chromosome   | Start     | End       |
  | (category)   | (int32)   | (int32)   |
  |--------------+-----------+-----------|
  | chr1         | 100079649 | 100079674 |
  | chr1         | 212609534 | 212609559 |
  | chr1         | 223587418 | 223587443 |
  | chr1         | 202450161 | 202450186 |
  | ...          | ...       | ...       |
  | chrY         | 11942770  | 11942795  |
  | chrY         | 8316773   | 8316798   |
  | chrY         | 7463444   | 7463469   |
  | chrY         | 7405376   | 7405401   |
  +--------------+-----------+-----------+
  Unstranded PyRanges object has 10,000 rows and 3 columns from 24 chromosomes.
  For printing, the PyRanges was sorted on Chromosome.


It is possible to make PyRanges objects out of basic Python datatypes:

  >>> gr = pr.PyRanges(chromosomes="chr1", strands="+", starts=[0, 1, 2], ends=(3, 4, 5))
  >>> print(gr)
  +--------------+-----------+-----------+--------------+
  | Chromosome   |     Start |       End | Strand       |
  | (category)   |   (int32) |   (int32) | (category)   |
  |--------------+-----------+-----------+--------------|
  | chr1         |         0 |         3 | +            |
  | chr1         |         1 |         4 | +            |
  | chr1         |         2 |         5 | +            |
  +--------------+-----------+-----------+--------------+
  Stranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

  >>> gr = pr.PyRanges(chromosomes="chr1 chr2 chr3".split(), strands="+ - +".split(), starts=[0, 1, 2], ends=(3, 4, 5))
  >>> print(gr)
  +--------------+-----------+-----------+--------------+
  | Chromosome   |     Start |       End | Strand       |
  | (category)   |   (int32) |   (int32) | (category)   |
  |--------------+-----------+-----------+--------------|
  | chr1         |         0 |         3 | +            |
  | chr2         |         1 |         4 | -            |
  | chr3         |         2 |         5 | +            |
  +--------------+-----------+-----------+--------------+
  Stranded PyRanges object has 3 rows and 4 columns from 3 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
	
Using read_bed, read_gtf, read_gff3 or read_bam
-----------------------------------------------


The pyranges library can create PyRanges from gff3 common file formats, namely gtf/gff, gff3, bed and bam ^.

  >>> ensembl_path = pr.get_example_path("ensembl.gtf")
  >>> gr = pr.read_gtf(ensembl_path)
  >>> print(gr)
  +--------------+------------+--------------+-----------+-----------+------------+-------+
  | Chromosome   | Source     | Feature      | Start     | End       | Score      | +20   |
  | (category)   | (object)   | (category)   | (int32)   | (int32)   | (object)   | ...   |
  |--------------+------------+--------------+-----------+-----------+------------+-------|
  | 1            | havana     | gene         | 11868     | 14409     | .          | ...   |
  | 1            | havana     | transcript   | 11868     | 14409     | .          | ...   |
  | 1            | havana     | exon         | 11868     | 12227     | .          | ...   |
  | 1            | havana     | exon         | 12612     | 12721     | .          | ...   |
  | ...          | ...        | ...          | ...       | ...       | ...        | ...   |
  | 1            | ensembl    | transcript   | 120724    | 133723    | .          | ...   |
  | 1            | ensembl    | exon         | 133373    | 133723    | .          | ...   |
  | 1            | ensembl    | exon         | 129054    | 129223    | .          | ...   |
  | 1            | ensembl    | exon         | 120873    | 120932    | .          | ...   |
  +--------------+------------+--------------+-----------+-----------+------------+-------+
  Stranded PyRanges object has 95 rows and 26 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  20 hidden columns: Strand, Frame, gene_id, gene_version, gene_name, gene_source, gene_biotype, ... (+ 13 more.)


To read bam files the optional bamread-library must be installed. Use::

    conda install -c bioconda bamread
 
or::
    
    pip install bamread 

to install it
    
    
read_bam takes the arguments ``sparse``, ``mapq``, ``required_flag``, ``filter_flag``, which have the default values True, 0, 0 and 1540, respectively. With sparse True, only the columns ``['Chromosome', 'Start', 'End', 'Strand', 'Flag']`` are fetched. Setting sparse to False additionally gives you the columns ``['QueryStart', 'QueryEnd', 'Name', 'Cigar', 'Quality']``, but is more time and memory-consuming.
All the reader functions also take the flag ``as_df``


Using from_dict
---------------

  >>> f1 = pr.data.f1()
  >>> d = f1.to_example(n=10)
  >>> print(d)
  {'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [3, 8, 5], 'End': [6, 9, 7], 'Name': ['interval1', 'interval3', 'interval2'], 'Score': [0, 0, 0], 'Strand': ['+', '+', '-']}
	
  >>> print(pr.from_dict(d))
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   |     Start |       End | Name       |     Score | Strand       |
  | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chr1         |         3 |         6 | interval1  |         0 | +            |
  | chr1         |         8 |         9 | interval3  |         0 | +            |
  | chr1         |         5 |         7 | interval2  |         0 | -            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.


Writing PyRanges to disk
~~~~~~~~~~~~~~~~~~~~~~~~


The PyRanges can be written to several formats, namely csv, gtf, gff3 and bigwig.
If no path-argument is given, the string representation of the data is returned. (It may potentially be very large.) If a path is given, the return value is the object itself. This way the write methods can easily be inserted in method call chains.

  >>> import pyranges as pr
  >>> gr = pr.data.chipseq()
  >>> gr.to_gtf("chipseq.gtf")
  # file chipseq.gtf has been created 



The to_csv method takes the arguments header and sep.

  >>> print(gr.drop(['Label', 'Tag']).head().to_csv(sep="\t", header=False))
  chr1	212609534	212609559	U0	0	+
  chr1	169887529	169887554	U0	0	+
  chr1	216711011	216711036	U0	0	+
  chr1	144227079	144227104	U0	0	+
  chr1	148177825	148177850	U0	0	+
  chr1	113486652	113486677	U0	0	+
  chr1	27024083	27024108	U0	0	+
  chr1	37865066	37865091	U0	0	+

All to-methods except to_bigwig takes an argument chain which can be set to True if you want the method to return the PyRanges it wrote. It is useful for storing the intermediate results of long call chains.::

	pr.data().f1().to_csv("bla", chain=True).merge()...
	
	
	
The pyranges library can also create bigwigs, but it needs the library pybigwig which is not installed by default. 
Use:: 
	
	conda install -c bioconda pybigwig
	
or::

	pip install pybigwig
	

to install it.

The bigwig writer needs to know the chromosome sizes. 
You can fetch these using the pyranges database functions, a pyranges add-on that can be install with:

.. code-block:: bash

	pip install pyranges_db

.. doctest::

  >>> gr.to_bigwig("chipseq.bw", chromsizes)
  # file chipseq.bw has been created 



To create a bigwig from an arbitrary value column, use the value_col argument.
If you want to write one bigwig for each strand, you need to do it manually.

  >>> gr["+"].to_bigwig("chipseq_plus.bw", chromsizes)
  >>> gr["-"].to_bigwig("chipseq_minus.bw", chromsizes)

to_bigwig also takes a flag ``divide_by`` which takes another PyRanges. Using divide_by creates a log2-normalized bigwig.




Inspecting PyRanges
~~~~~~~~~~~~~~~~~~~


The PyRanges method print provides an overview of its data:


  >>> import pyranges as pr
  >>> gr = pr.data.chipseq()
  >>> gr.print()
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   | Start     | End       | Name       | Score     | Strand       |
  | (category)   | (int32)   | (int32)   | (object)   | (int64)   | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chr1         | 212609534 | 212609559 | U0         | 0         | +            |
  | chr1         | 169887529 | 169887554 | U0         | 0         | +            |
  | chr1         | 216711011 | 216711036 | U0         | 0         | +            |
  | chr1         | 144227079 | 144227104 | U0         | 0         | +            |
  | ...          | ...       | ...       | ...        | ...       | ...          |
  | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
  | chrY         | 13517892  | 13517917  | U0         | 0         | -            |
  | chrY         | 8010951   | 8010976   | U0         | 0         | -            |
  | chrY         | 7405376   | 7405401   | U0         | 0         | -            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 10,000 rows and 6 columns from 24 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

The same method is invoked under the hood anytime we request a string representation:

  >>> print(str(gr))
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   | Start     | End       | Name       | Score     | Strand       |
  | (category)   | (int32)   | (int32)   | (object)   | (int64)   | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chr1         | 212609534 | 212609559 | U0         | 0         | +            |
  | chr1         | 169887529 | 169887554 | U0         | 0         | +            |
  | chr1         | 216711011 | 216711036 | U0         | 0         | +            |
  | chr1         | 144227079 | 144227104 | U0         | 0         | +            |
  | ...          | ...       | ...       | ...        | ...       | ...          |
  | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
  | chrY         | 13517892  | 13517917  | U0         | 0         | -            |
  | chrY         | 8010951   | 8010976   | U0         | 0         | -            |
  | chrY         | 7405376   | 7405401   | U0         | 0         | -            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 10,000 rows and 6 columns from 24 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

As explained in the tutorial, PyRanges objects consist of collections of DataFrames, organized per chromosome (and strand, if Stranded). When printed, they are displayed as a continuous table, ordered by Chromosome (and strand). 

The window width affects the output of print: columns that do not fit are hidden. When this happens, a message is printed after the table:

  >>> gr.new_col = 'value'
  >>> gr.another_col = 99
  >>> gr.print()
  +--------------+-----------+-----------+------------+-----------+-------+
  | Chromosome   | Start     | End       | Name       | Score     | +3    |
  | (category)   | (int32)   | (int32)   | (object)   | (int64)   | ...   |
  |--------------+-----------+-----------+------------+-----------+-------|
  | chr1         | 212609534 | 212609559 | U0         | 0         | ...   |
  | chr1         | 169887529 | 169887554 | U0         | 0         | ...   |
  | chr1         | 216711011 | 216711036 | U0         | 0         | ...   |
  | chr1         | 144227079 | 144227104 | U0         | 0         | ...   |
  | ...          | ...       | ...       | ...        | ...       | ...   |
  | chrY         | 15224235  | 15224260  | U0         | 0         | ...   |
  | chrY         | 13517892  | 13517917  | U0         | 0         | ...   |
  | chrY         | 8010951   | 8010976   | U0         | 0         | ...   |
  | chrY         | 7405376   | 7405401   | U0         | 0         | ...   |
  +--------------+-----------+-----------+------------+-----------+-------+
  Stranded PyRanges object has 10,000 rows and 8 columns from 24 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  3 hidden columns: Strand, new_col, another_col

Only a limited number of rows are displayed, which are taken from the top and bottom of the table. This is 8 by default, and can be redefined through the first argument of print, named n:

  >>> gr.print(2)
  +--------------+-----------+-----------+------------+-----------+-------+
  | Chromosome   | Start     | End       | Name       | Score     | +3    |
  | (category)   | (int32)   | (int32)   | (object)   | (int64)   | ...   |
  |--------------+-----------+-----------+------------+-----------+-------|
  | chr1         | 212609534 | 212609559 | U0         | 0         | ...   |
  | ...          | ...       | ...       | ...        | ...       | ...   |
  | chrY         | 7405376   | 7405401   | U0         | 0         | ...   |
  +--------------+-----------+-----------+------------+-----------+-------+
  Stranded PyRanges object has 10,000 rows and 8 columns from 24 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  3 hidden columns: Strand, new_col, another_col

  >>> gr.print(n=20)
  +--------------+-----------+-----------+------------+-----------+-------+
  | Chromosome   | Start     | End       | Name       | Score     | +3    |
  | (category)   | (int32)   | (int32)   | (object)   | (int64)   | ...   |
  |--------------+-----------+-----------+------------+-----------+-------|
  | chr1         | 212609534 | 212609559 | U0         | 0         | ...   |
  | chr1         | 169887529 | 169887554 | U0         | 0         | ...   |
  | chr1         | 216711011 | 216711036 | U0         | 0         | ...   |
  | chr1         | 144227079 | 144227104 | U0         | 0         | ...   |
  | chr1         | 148177825 | 148177850 | U0         | 0         | ...   |
  | chr1         | 113486652 | 113486677 | U0         | 0         | ...   |
  | chr1         | 27024083  | 27024108  | U0         | 0         | ...   |
  | chr1         | 37865066  | 37865091  | U0         | 0         | ...   |
  | chr1         | 47488200  | 47488225  | U0         | 0         | ...   |
  | chr1         | 197075093 | 197075118 | U0         | 0         | ...   |
  | ...          | ...       | ...       | ...        | ...       | ...   |
  | chrY         | 21707662  | 21707687  | U0         | 0         | ...   |
  | chrY         | 7761026   | 7761051   | U0         | 0         | ...   |
  | chrY         | 22210637  | 22210662  | U0         | 0         | ...   |
  | chrY         | 14774053  | 14774078  | U0         | 0         | ...   |
  | chrY         | 16495497  | 16495522  | U0         | 0         | ...   |
  | chrY         | 7046809   | 7046834   | U0         | 0         | ...   |
  | chrY         | 15224235  | 15224260  | U0         | 0         | ...   |
  | chrY         | 13517892  | 13517917  | U0         | 0         | ...   |
  | chrY         | 8010951   | 8010976   | U0         | 0         | ...   |
  | chrY         | 7405376   | 7405401   | U0         | 0         | ...   |
  +--------------+-----------+-----------+------------+-----------+-------+
  Stranded PyRanges object has 10,000 rows and 8 columns from 24 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  3 hidden columns: Strand, new_col, another_col

Argument formatting allows to fine-tune appearance. It takes a dictionary with any column name as key, and a string as value which follows the python format syntax:

  >>> gr.print(formatting={
  ...     'Score':'{:.2f}',
  ... 	    'End':'{:e}',
  ... 	    'Start':'{:,}',
  ... 	    'Name':'name={}',
  ... 	    })
  +--------------+-------------+--------------+------------+-----------+-------+
  | Chromosome   | Start       | End          | Name       | Score     | +3    |
  | (category)   | (int32)     | (int32)      | (object)   | (int64)   | ...   |
  |--------------+-------------+--------------+------------+-----------+-------|
  | chr1         | 212,609,534 | 2.126096e+08 | name=U0    | 0.00      | ...   |
  | chr1         | 169,887,529 | 1.698876e+08 | name=U0    | 0.00      | ...   |
  | chr1         | 216,711,011 | 2.167110e+08 | name=U0    | 0.00      | ...   |
  | chr1         | 144,227,079 | 1.442271e+08 | name=U0    | 0.00      | ...   |
  | ...          | ...         | ...          | ...        | ...       | ...   |
  | chrY         | 15,224,235  | 1.522426e+07 | name=U0    | 0.00      | ...   |
  | chrY         | 13,517,892  | 1.351792e+07 | name=U0    | 0.00      | ...   |
  | chrY         | 8,010,951   | 8.010976e+06 | name=U0    | 0.00      | ...   |
  | chrY         | 7,405,376   | 7.405401e+06 | name=U0    | 0.00      | ...   |
  +--------------+-------------+--------------+------------+-----------+-------+
  Stranded PyRanges object has 10,000 rows and 8 columns from 24 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  3 hidden columns: Strand, new_col, another_col


PyRanges columns are pandas Series, and they may be of different data types. The types  are shown in the header shown with print (see above). To see them all, use property dtypes:

  >>> gr.dtypes
  Chromosome     category
  Start             int32
  End               int32
  Name             object
  Score             int64
  Strand         category
  new_col          object
  another_col       int64
  dtype: object

If you want to inspect more information from a PyRanges object, remember that you can always transform it into a pandas DataFrame, which gives access to all its methods. For example, you may employ pandas info and describe:

  >>> gr.df.info()
  <class 'pandas.core.frame.DataFrame'>
  RangeIndex: 10000 entries, 0 to 9999
  Data columns (total 8 columns):
      Column       Non-Null Count  Dtype
   ---  ------       --------------  -----
   0   Chromosome   10000 non-null  category
   1   Start        10000 non-null  int32
   2   End          10000 non-null  int32
   3   Name         10000 non-null  object
   4   Score        10000 non-null  int64
   5   Strand       10000 non-null  category
   6   new_col      10000 non-null  object
   7   another_col  10000 non-null  int64
   dtypes: category(2), int32(2), int64(2), object(2)
   memory usage: 411.1+ KB

  >>> gr.df.describe()
        Start           End    Score  another_col
  count  1.000000e+04  1.000000e+04  10000.0      10000.0
  mean   8.087570e+07  8.087573e+07      0.0         99.0
  std    5.572825e+07  5.572825e+07      0.0          0.0
  min    1.361100e+04  1.363600e+04      0.0         99.0
  25%    3.550257e+07  3.550260e+07      0.0         99.0
  50%    7.030672e+07  7.030674e+07      0.0         99.0
  75%    1.167902e+08  1.167902e+08      0.0         99.0
  max    2.471349e+08  2.471349e+08      0.0         99.0


Accessing data
~~~~~~~~~~~~~~

Selecting rows
--------------

As seen in the tutorial, PyRanges provides various ways to select a subset of rows. All of these methods return a (smaller) copy of the original object.

One way is to index **by genomic region**, which may take any of the following syntaxes:

* chromosome
* chromosome, position slice 
* chromosome, strand, position slice

Here's one example for each:

.. code-block:: python

  >>> import pyranges as pr
  >>> gr = pr.data.chipseq()
  >>> gr['chrX']
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   | Start     | End       | Name       | Score     | Strand       |
  | (category)   | (int32)   | (int32)   | (object)   | (int64)   | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chrX         | 13843759  | 13843784  | U0         | 0         | +            |
  | chrX         | 114673546 | 114673571 | U0         | 0         | +            |
  | chrX         | 131816774 | 131816799 | U0         | 0         | +            |
  | chrX         | 45504745  | 45504770  | U0         | 0         | +            |
  | ...          | ...       | ...       | ...        | ...       | ...          |
  | chrX         | 146694149 | 146694174 | U0         | 0         | -            |
  | chrX         | 5044527   | 5044552   | U0         | 0         | -            |
  | chrX         | 15281263  | 15281288  | U0         | 0         | -            |
  | chrX         | 120273723 | 120273748 | U0         | 0         | -            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 282 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

  >>> gr['chr1', 1000000:3000000]
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   |     Start |       End | Name       |     Score | Strand       |
  | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chr1         |   1541598 |   1541623 | U0         |         0 | +            |
  | chr1         |   1599121 |   1599146 | U0         |         0 | +            |
  | chr1         |   1325303 |   1325328 | U0         |         0 | -            |
  | chr1         |   1820285 |   1820310 | U0         |         0 | -            |
  | chr1         |   2448322 |   2448347 | U0         |         0 | -            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 5 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

  >>> gr['chr1', '-', 1000000:3000000]
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   |     Start |       End | Name       |     Score | Strand       |
  | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chr1         |   1325303 |   1325328 | U0         |         0 | -            |
  | chr1         |   1820285 |   1820310 | U0         |         0 | -            |
  | chr1         |   2448322 |   2448347 | U0         |         0 | -            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

Simple forms of row selection are done through methods **head** and **tail**, which return the top or bottom N rows, respectively, where N is 8 by default:

  >>> gr.head()
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   |     Start |       End | Name       |     Score | Strand       |
  | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chr1         | 212609534 | 212609559 | U0         |         0 | +            |
  | chr1         | 169887529 | 169887554 | U0         |         0 | +            |
  | chr1         | 216711011 | 216711036 | U0         |         0 | +            |
  | chr1         | 144227079 | 144227104 | U0         |         0 | +            |
  | chr1         | 148177825 | 148177850 | U0         |         0 | +            |
  | chr1         | 113486652 | 113486677 | U0         |         0 | +            |
  | chr1         |  27024083 |  27024108 | U0         |         0 | +            |
  | chr1         |  37865066 |  37865091 | U0         |         0 | +            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 8 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

  >>> gr.tail()
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   |     Start |       End | Name       |     Score | Strand       |
  | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chrY         |  22210637 |  22210662 | U0         |         0 | -            |
  | chrY         |  14774053 |  14774078 | U0         |         0 | -            |
  | chrY         |  16495497 |  16495522 | U0         |         0 | -            |
  | chrY         |   7046809 |   7046834 | U0         |         0 | -            |
  | chrY         |  15224235 |  15224260 | U0         |         0 | -            |
  | chrY         |  13517892 |  13517917 | U0         |         0 | -            |
  | chrY         |   8010951 |   8010976 | U0         |         0 | -            |
  | chrY         |   7405376 |   7405401 | U0         |         0 | -            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 8 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.


The most important form of row selection is by **indexing with a boolean Series**. This is typically generated from a column through a comparison operator. Let's see it with some other example data:

  >>> gg = pr.data.chipseq()
  >>> gg.print(n=20)
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   |     Start |       End | Name       |     Score | Strand       |
  | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chr1         |      9939 |     10138 | H3K27me3   |         7 | +            |
  | chr1         |      9953 |     10152 | H3K27me3   |         5 | +            |
  | chr1         |     10024 |     10223 | H3K27me3   |         1 | +            |
  | chr1         |     10246 |     10445 | H3K27me3   |         4 | +            |
  | chr1         |    110246 |    110445 | H3K27me3   |         1 | +            |
  | chr1         |      9916 |     10115 | H3K27me3   |         5 | -            |
  | chr1         |      9951 |     10150 | H3K27me3   |         8 | -            |
  | chr1         |      9978 |     10177 | H3K27me3   |         7 | -            |
  | chr1         |     10001 |     10200 | H3K27me3   |         5 | -            |
  | chr1         |     10127 |     10326 | H3K27me3   |         1 | -            |
  | chr1         |     10241 |     10440 | H3K27me3   |         6 | -            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 11 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

Below, we produce a boolean Series:

  >>> gg.Score > 5
  1      True
  3     False
  6     False
  9     False
  10    False
  0     False
  2      True
  4      True
  5     False
  7     False
  8      True
  Name: Score, dtype: bool

And we use it to select the rows in which the column Score has a value greater than 5:

  >>> gg[gg.Score>5]
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   |     Start |       End | Name       |     Score | Strand       |
  | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chr1         |      9939 |     10138 | H3K27me3   |         7 | +            |
  | chr1         |      9951 |     10150 | H3K27me3   |         8 | -            |
  | chr1         |      9978 |     10177 | H3K27me3   |         7 | -            |
  | chr1         |     10241 |     10440 | H3K27me3   |         6 | -            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 4 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

As pandas users know, these logical operators can be employed with boolean Series:

* "&" =  element-wise AND operator
* "|" = element-wise OR operator
* "~" = NOT operator, inverts the values of the Series on its right

When using logical operators, make sure to parenthesize properly. 

Let's get the + intervals with Score 1 starting before 12,000 or ending after 100,000:

  >>> gg[
  ...    (gg.Score==1) &
  ...    (gg.Strand=='+') &
  ...    ((gg.Start<12000) | (gg.End>100000))
  ...    ]
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   |     Start |       End | Name       |     Score | Strand       |
  | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chr1         |     10024 |     10223 | H3K27me3   |         1 | +            |
  | chr1         |    110246 |    110445 | H3K27me3   |         1 | +            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

Let's invert the selection, i.e. taking all intervals that do not fit the above criteria:

  >>> gg[~(
  ...      (gg.Score==1) &
  ...      (gg.Strand=='+') &
  ...      ((gg.Start<12000) | (gg.End>100000))
  ...     )
  ...    ]
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   | Start     | End       | Name       | Score     | Strand       |
  | (category)   | (int32)   | (int32)   | (object)   | (int64)   | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chr1         | 9939      | 10138     | H3K27me3   | 7         | +            |
  | chr1         | 9953      | 10152     | H3K27me3   | 5         | +            |
  | chr1         | 10246     | 10445     | H3K27me3   | 4         | +            |
  | chr1         | 9916      | 10115     | H3K27me3   | 5         | -            |
  | ...          | ...       | ...       | ...        | ...       | ...          |
  | chr1         | 9978      | 10177     | H3K27me3   | 7         | -            |
  | chr1         | 10001     | 10200     | H3K27me3   | 5         | -            |
  | chr1         | 10127     | 10326     | H3K27me3   | 1         | -            |
  | chr1         | 10241     | 10440     | H3K27me3   | 6         | -            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 9 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.


Another way to select rows is **the subset method**, for which you provide a function which is then applied to each DataFrame of the collection, and which must return a boolean Series. Typically, you define a lambda function on-the-fly:


  
  >>> # the following is equivalent to
  >>> gg[gg.Score.isin([2,4,6]]
  >>> gg.subset(lambda x:x.Score.isin([2,4,6]))
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   |     Start |       End | Name       |     Score | Strand       |
  | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chr1         |     10246 |     10445 | H3K27me3   |         4 | +            |
  | chr1         |     10241 |     10440 | H3K27me3   |         6 | -            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

The method *subset* is suited for complex pandas operations, and it is also useful in method call chains. 

Lastly, a fairly specific form of row selection is **drop_duplicate_positions**, which gets rid of interval with the same coordinates:

  >>> d = {"Chromosome": [1, 1, 1, 2, 2], 
  ...      "Start": [1, 1, 2, 1, 8], 
  ...      "End": [4, 4, 9, 4, 12], 
  ...      "Strand": ["+", "+", "+", "+","-"], 
  ...      "ID": ["a", "b", "c", "d", "e"]}
  >>> p = pr.from_dict(d)
  >>> p
  +--------------+-----------+-----------+--------------+------------+
  |   Chromosome |     Start |       End | Strand       | ID         |
  |   (category) |   (int32) |   (int32) | (category)   | (object)   |
  |--------------+-----------+-----------+--------------+------------|
  |            1 |         1 |         4 | +            | a          |
  |            1 |         1 |         4 | +            | b          |
  |            1 |         2 |         9 | +            | c          |
  |            2 |         1 |         4 | +            | d          |
  |            2 |         8 |        12 | -            | e          |
  +--------------+-----------+-----------+--------------+------------+
  Stranded PyRanges object has 5 rows and 5 columns from 2 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
	
	
  >>> p.drop_duplicate_positions()
  +--------------+-----------+-----------+--------------+------------+
  |   Chromosome |     Start |       End | Strand       | ID         |
  |   (category) |   (int32) |   (int32) | (category)   | (object)   |
  |--------------+-----------+-----------+--------------+------------|
  |            1 |         1 |         4 | +            | a          |
  |            1 |         2 |         9 | +            | c          |
  |            2 |         1 |         4 | +            | d          |
  |            2 |         8 |        12 | -            | e          |
  +--------------+-----------+-----------+--------------+------------+
  Stranded PyRanges object has 4 rows and 5 columns from 2 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

Normally, the first instance of duplicated intervals is retained. Through argument keep=False, you can decide to remove them all:

  >>> p.drop_duplicate_positions(keep=False)
  +--------------+-----------+-----------+--------------+------------+
  |   Chromosome |     Start |       End | Strand       | ID         |
  |   (category) |   (int32) |   (int32) | (category)   | (object)   |
  |--------------+-----------+-----------+--------------+------------|
  |            1 |         2 |         9 | +            | c          |
  |            2 |         1 |         4 | +            | d          |
  |            2 |         8 |        12 | -            | e          |
  +--------------+-----------+-----------+--------------+------------+
  Stranded PyRanges object has 3 rows and 5 columns from 2 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

Selecting columns
-----------------


As previously seen, single PyRanges column (which are pandas Series) can be extracted through the dot notation:


  >>> gr = pr.data.chipseq()
  >>> gr.Chromosome
  18      chr1
  70      chr1
  129     chr1
  170     chr1
  196     chr1
  	...
  3023    chrY
  3131    chrY
  3816    chrY
  3897    chrY
  9570    chrY
  Name: Chromosome, Length: 10000, dtype: category
  Categories (24, object): ['chr1', 'chr10', 'chr11', 'chr12', ..., 'chr8', 'chr9', 'chrX', 'chrY']

The same syntax can be used for the core PyRanges columns (Chromosome, Strand, Start, End) or for metadata columns:

  >>> gr.Name
  18      U0
  70      U0
  129     U0
  170     U0
  196     U0
  	...
  3023    U0
  3131    U0
  3816    U0
  3897    U0
  9570    U0
  Name: Name, Length: 10000, dtype: object

This syntax is analogous to pandas Dataframes. Note that, however, the bracket column selection in pandas does not work in the same way in PyRanges:

  >>> df=gr.df
  >>> df['End']
  0       212609559
  1       169887554
  2       216711036
  3       144227104
  4       148177850
  	  ...
  9995      7046834
  9996     15224260
  9997     13517917
  9998      8010976
  9999      7405401
  Name: End, Length: 10000, dtype: int32

  >>> gr['End']
  Empty PyRanges

Because the last expression is evaluated as a genomic region, i.e. a form of row selection: it is searching for intervals on a Chromosome named "End", and finds none. Indeed, this fetches intervals on the chrY:

  >>> gr['chrY']
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   | Start     | End       | Name       | Score     | Strand       |
  | (category)   | (int32)   | (int32)   | (object)   | (int64)   | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chrY         | 12930373  | 12930398  | U0         | 0         | +            |
  | chrY         | 15548022  | 15548047  | U0         | 0         | +            |
  | chrY         | 7194340   | 7194365   | U0         | 0         | +            |
  | chrY         | 21559181  | 21559206  | U0         | 0         | +            |
  | ...          | ...       | ...       | ...        | ...       | ...          |
  | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
  | chrY         | 13517892  | 13517917  | U0         | 0         | -            |
  | chrY         | 8010951   | 8010976   | U0         | 0         | -            |
  | chrY         | 7405376   | 7405401   | U0         | 0         | -            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 23 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

You can provide a list of column names in the bracket notation to select those columns. Pyranges will still return a PyRanges object, therefore retaining the core columns regardless of whether they were selected or not:

  >>> gr[ ['Name'] ]
  +--------------+-----------+-----------+------------+--------------+
  | Chromosome   | Start     | End       | Name       | Strand       |
  | (category)   | (int32)   | (int32)   | (object)   | (category)   |
  |--------------+-----------+-----------+------------+--------------|
  | chr1         | 212609534 | 212609559 | U0         | +            |
  | chr1         | 169887529 | 169887554 | U0         | +            |
  | chr1         | 216711011 | 216711036 | U0         | +            |
  | chr1         | 144227079 | 144227104 | U0         | +            |
  | ...          | ...       | ...       | ...        | ...          |
  | chrY         | 15224235  | 15224260  | U0         | -            |
  | chrY         | 13517892  | 13517917  | U0         | -            |
  | chrY         | 8010951   | 8010976   | U0         | -            |
  | chrY         | 7405376   | 7405401   | U0         | -            |
  +--------------+-----------+-----------+------------+--------------+
  Stranded PyRanges object has 10,000 rows and 5 columns from 24 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

This is convenient to reduce genome annotation that consists of many columns:

  >>> ensembl_path = pr.get_example_path("ensembl.gtf")
  >>> ge = pr.read_gtf(ensembl_path)
  >>> ge
  +--------------+------------+--------------+-----------+-----------+------------+-------+
  | Chromosome   | Source     | Feature      | Start     | End       | Score      | +20   |
  | (category)   | (object)   | (category)   | (int32)   | (int32)   | (object)   | ...   |
  |--------------+------------+--------------+-----------+-----------+------------+-------|
  | 1            | havana     | gene         | 11868     | 14409     | .          | ...   |
  | 1            | havana     | transcript   | 11868     | 14409     | .          | ...   |
  | 1            | havana     | exon         | 11868     | 12227     | .          | ...   |
  | 1            | havana     | exon         | 12612     | 12721     | .          | ...   |
  | ...          | ...        | ...          | ...       | ...       | ...        | ...   |
  | 1            | ensembl    | transcript   | 120724    | 133723    | .          | ...   |
  | 1            | ensembl    | exon         | 133373    | 133723    | .          | ...   |
  | 1            | ensembl    | exon         | 129054    | 129223    | .          | ...   |
  | 1            | ensembl    | exon         | 120873    | 120932    | .          | ...   |
  +--------------+------------+--------------+-----------+-----------+------------+-------+
  Stranded PyRanges object has 95 rows and 26 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
	20 hidden columns: Strand, Frame, gene_id, gene_version, gene_name, gene_source, gene_biotype, ... (+ 13 more.)

  >>> ge[ ['gene_id', 'gene_name'] ]
  +--------------+-----------+-----------+--------------+-----------------+-------------+
  | Chromosome   | Start     | End       | Strand       | gene_id         | gene_name   |
  | (category)   | (int32)   | (int32)   | (category)   | (object)        | (object)    |
  |--------------+-----------+-----------+--------------+-----------------+-------------|
  | 1            | 11868     | 14409     | +            | ENSG00000223972 | DDX11L1     |
  | 1            | 11868     | 14409     | +            | ENSG00000223972 | DDX11L1     |
  | 1            | 11868     | 12227     | +            | ENSG00000223972 | DDX11L1     |
  | 1            | 12612     | 12721     | +            | ENSG00000223972 | DDX11L1     |
  | ...          | ...       | ...       | ...          | ...             | ...         |
  | 1            | 120724    | 133723    | -            | ENSG00000238009 | AL627309.1  |
  | 1            | 133373    | 133723    | -            | ENSG00000238009 | AL627309.1  |
  | 1            | 129054    | 129223    | -            | ENSG00000238009 | AL627309.1  |
  | 1            | 120873    | 120932    | -            | ENSG00000238009 | AL627309.1  |
  +--------------+-----------+-----------+--------------+-----------------+-------------+
  Stranded PyRanges object has 95 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

The **drop method** is an alternative way of column selection wherein we specify what we want to remove, rather than what to keep:


  >>> gr.print()
  >>> gr.drop(['Name']).print()
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   | Start     | End       | Name       | Score     | Strand       |
  | (category)   | (int32)   | (int32)   | (object)   | (int64)   | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chr1         | 212609534 | 212609559 | U0         | 0         | +            |
  | chr1         | 169887529 | 169887554 | U0         | 0         | +            |
  | chr1         | 216711011 | 216711036 | U0         | 0         | +            |
  | chr1         | 144227079 | 144227104 | U0         | 0         | +            |
  | ...          | ...       | ...       | ...        | ...       | ...          |
  | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
  | chrY         | 13517892  | 13517917  | U0         | 0         | -            |
  | chrY         | 8010951   | 8010976   | U0         | 0         | -            |
  | chrY         | 7405376   | 7405401   | U0         | 0         | -            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 10,000 rows and 6 columns from 24 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  +--------------+-----------+-----------+-----------+--------------+
  | Chromosome   | Start     | End       | Score     | Strand       |
  | (category)   | (int32)   | (int32)   | (int64)   | (category)   |
  |--------------+-----------+-----------+-----------+--------------|
  | chr1         | 212609534 | 212609559 | 0         | +            |
  | chr1         | 169887529 | 169887554 | 0         | +            |
  | chr1         | 216711011 | 216711036 | 0         | +            |
  | chr1         | 144227079 | 144227104 | 0         | +            |
  | ...          | ...       | ...       | ...       | ...          |
  | chrY         | 15224235  | 15224260  | 0         | -            |
  | chrY         | 13517892  | 13517917  | 0         | -            |
  | chrY         | 8010951   | 8010976   | 0         | -            |
  | chrY         | 7405376   | 7405401   | 0         | -            |
  +--------------+-----------+-----------+-----------+--------------+
  Stranded PyRanges object has 10,000 rows and 5 columns from 24 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

Without arguments, drop will get rid of all non-core columns:

  >>> gr.drop()
  +--------------+-----------+-----------+--------------+
  | Chromosome   | Start     | End       | Strand       |
  | (category)   | (int32)   | (int32)   | (category)   |
  |--------------+-----------+-----------+--------------|
  | chr1         | 212609534 | 212609559 | +            |
  | chr1         | 169887529 | 169887554 | +            |
  | chr1         | 216711011 | 216711036 | +            |
  | chr1         | 144227079 | 144227104 | +            |
  | ...          | ...       | ...       | ...          |
  | chrY         | 15224235  | 15224260  | -            |
  | chrY         | 13517892  | 13517917  | -            |
  | chrY         | 8010951   | 8010976   | -            |
  | chrY         | 7405376   | 7405401   | -            |
  +--------------+-----------+-----------+--------------+
  Stranded PyRanges object has 10,000 rows and 4 columns from 24 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.


If you want to obtain a DataFrame with certain columns rather than a PyRanges object, get a DataFrame copy through the df property, then perform pandas-style column selection. Obviously, in this case core columns are returned only if explicitly selected:

  >>> gr.df [ ['Name', 'Start'] ]
       Name      Start
  0      U0  212609534
  1      U0  169887529
  2      U0  216711011
  3      U0  144227079
  4      U0  148177825
  ...   ...        ...
  9995   U0    7046809
  9996   U0   15224235
  9997   U0   13517892
  9998   U0    8010951
  9999   U0    7405376
  
  [10000 rows x 2 columns]



Obtaining sequences
-------------------


A common operation is to fetch the sequences corresponding to the intervals represented in the PyRanges object. Function ``get_sequence`` takes as input a PyRanges object and the path to a fasta file, and returns a Series containing sequences, in the same order as the intervals. It requires package pyfaidx (install with pip install pyfaidx).

In the tutorial, we saw its usage with a real genome. Let's make a toy example here:

  >>> with open('minigenome.fa', 'w') as fw:
  ...     fw.write('>chrZ\n')
  ...     fw.write('AAAGGGCCCTTTAAAGGGCCCTTTAAAGGGCCCTTT\n')

  >>> sg = pr.from_dict({"Chromosome": ["chrZ", "chrZ", "chrZ", "chrZ"],
  ... 	           "Start": [0, 5, 10, 10], "End": [3, 8, 20, 20],
  ... 	           "name":["a", "a", "b", "c"],
  ... 	           "Strand":["+", "+", "+", "-"] })
  
  >>> sg 
  +--------------+-----------+-----------+------------+--------------+
  | Chromosome   |     Start |       End | name       | Strand       |
  | (category)   |   (int32) |   (int32) | (object)   | (category)   |
  |--------------+-----------+-----------+------------+--------------|
  | chrZ         |         0 |         3 | a          | +            |
  | chrZ         |         5 |         8 | a          | +            |
  | chrZ         |        10 |        20 | b          | +            |
  | chrZ         |        10 |        20 | c          | -            |
  +--------------+-----------+-----------+------------+--------------+
  Stranded PyRanges object has 4 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

Note the genome sequence in the code above. Let's run ``get_sequences`` to obtain the portions corresponding to our intervals:


  >>> pr.get_sequence(sg, 'minigenome.fa')
  0           AAA
  1           GCC
  2    TTAAAGGGCC
  3    GGCCCTTTAA
  dtype: object

Note that the last two intervals have identical coordinates but are on opposite strands. Function ``get_sequence`` returns the reverse complement for intervals on the negative strand.

Since the returned Series has the same length as the PyRanges object, we can assign it to a new column:


  >>> sg.Sequence = pr.get_sequence(sg, 'minigenome.fa')
  >>> sg
  +--------------+-----------+-----------+------------+--------------+------------+
  | Chromosome   |     Start |       End | name       | Strand       | Sequence   |
  | (category)   |   (int32) |   (int32) | (object)   | (category)   | (object)   |
  |--------------+-----------+-----------+------------+--------------+------------|
  | chrZ         |         0 |         3 | a          | +            | AAA        |
  | chrZ         |         5 |         8 | a          | +            | GCC        |
  | chrZ         |        10 |        20 | b          | +            | TTAAAGGGCC |
  | chrZ         |        10 |        20 | c          | -            | GGCCCTTTAA |
  +--------------+-----------+-----------+------------+--------------+------------+
  Stranded PyRanges object has 4 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

This allows us to filter by sequence, using pandas string methods. For example, let's get those that start with G:



  >>> sg[sg.Sequence.str.startswith('G')]
  +--------------+-----------+-----------+------------+--------------+------------+
  | Chromosome   |     Start |       End | name       | Strand       | Sequence   |
  | (category)   |   (int32) |   (int32) | (object)   | (category)   | (object)   |
  |--------------+-----------+-----------+------------+--------------+------------|
  | chrZ         |         5 |         8 | a          | +            | GCC        |
  | chrZ         |        10 |        20 | c          | -            | GGCCCTTTAA |
  +--------------+-----------+-----------+------------+--------------+------------+
  Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

Let's get those which contain a CC and AA dinucleotides separated by 1-3 nucleotides:



  >>> sg[sg.Sequence.str.contains(r'CC.{1,3}AA', regex=True)]



Function ``get_sequence`` will treat each interval independently. Often, you want to get the sequence of an mRNA, i.e. concatenating exons. Function get_transcript_sequence serves this purpose, and employs argument group_by to group the exons into mRNAs:


  >>> pr.get_transcript_sequence(sg, group_by='name', path='minigenome.fa')
    name    Sequence
  0    a      AAAGCC
  1    b  TTAAAGGGCC
  2    c  GGCCCTTTAA

Note that this returns a pandas DataFrame with a row per exon group: its shape is different from the original PyRanges.



Operating with data
~~~~~~~~~~~~~~~~~~~


In this section, we give an overview of methods to modify the data in PyRanges.
Changing row order
Methods sort allows to sort intervals, i.e. altering the order of rows in the PyRanges object. When run without arguments, orders interval by increasing Start. Commonly, genomic annotation files are sorted in this way.


  >>> sg = pr.from_dict({"Chromosome": ["chrA", "chrA", "chrB", "chrB", "chrB"],
  ... 	           "Start": [55, 20, 65, 35, 75], 
  ... 	           "End": [88, 30, 75, 45, 85],
  ... 	           "name":["a", "a", "b", "c", "c"],
  ... 	           "Strand":["+", "+", "+", "-", "-"] })
  >>> sg
  +--------------+-----------+-----------+------------+--------------+
  | Chromosome   |     Start |       End | name       | Strand       |
  | (category)   |   (int32) |   (int32) | (object)   | (category)   |
  |--------------+-----------+-----------+------------+--------------|
  | chrA         |        55 |        88 | a          | +            |
  | chrA         |        20 |        30 | a          | +            |
  | chrB         |        65 |        75 | b          | +            |
  | chrB         |        35 |        45 | c          | -            |
  | chrB         |        75 |        85 | c          | -            |
  +--------------+-----------+-----------+------------+--------------+
  Stranded PyRanges object has 5 rows and 5 columns from 2 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

  >>> sg.sort()
  +--------------+-----------+-----------+------------+--------------+
  | Chromosome   |     Start |       End | name       | Strand       |
  | (category)   |   (int32) |   (int32) | (object)   | (category)   |
  |--------------+-----------+-----------+------------+--------------|
  | chrA         |        20 |        30 | a          | +            |
  | chrA         |        55 |        88 | a          | +            |
  | chrB         |        65 |        75 | b          | +            |
  | chrB         |        35 |        45 | c          | -            |
  | chrB         |        75 |        85 | c          | -            |
  +--------------+-----------+-----------+------------+--------------+
  Stranded PyRanges object has 5 rows and 5 columns from 2 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.


Remember that **sorting is performed separately for each internal table**: intervals on different chromosome/strands won't ever cross each other. To have all intervals sorted, work with a DataFrame object instead.

For intervals on the negative strand, it may be convenient to sort in the opposite order, since for them the leftmost exon is actually the last one in the mRNA. Instead of having to split the PyRanges object for this task, you may run sort with special argument "5", which will sort intervals in 5' to 3' order:


  >>> sg.sort('5')
  +--------------+-----------+-----------+------------+--------------+
  | Chromosome   |     Start |       End | name       | Strand       |
  | (category)   |   (int32) |   (int32) | (object)   | (category)   |
  |--------------+-----------+-----------+------------+--------------|
  | chrA         |        20 |        30 | a          | +            |
  | chrA         |        55 |        88 | a          | +            |
  | chrB         |        65 |        75 | b          | +            |
  | chrB         |        75 |        85 | c          | -            |
  | chrB         |        35 |        45 | c          | -            |
  +--------------+-----------+-----------+------------+--------------+
  Stranded PyRanges object has 5 rows and 5 columns from 2 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

Sorting may also take any column name, or a list of colum names, to sort rows by their value:

  >>> ag = pr.from_dict({"Chromosome": "chrX",
  ... 	           "Start": [55, 65, 20, 35, 75], 
  ... 	           "End": [88, 75, 30, 45, 85],
  ... 	           "Strand":["+", "+", "+", "+", "+"],
  ... 	           "col1":[1, 4, 4, 2, 2],
  ... 	            })
  >>> ag
  +--------------+-----------+-----------+--------------+-----------+
  | Chromosome   |     Start |       End | Strand       |      col1 |
  | (category)   |   (int32) |   (int32) | (category)   |   (int64) |
  |--------------+-----------+-----------+--------------+-----------|
  | chrX         |        55 |        88 | +            |         1 |
  | chrX         |        65 |        75 | +            |         4 |
  | chrX         |        20 |        30 | +            |         4 |
  | chrX         |        35 |        45 | +            |         2 |
  | chrX         |        75 |        85 | +            |         2 |
  +--------------+-----------+-----------+--------------+-----------+
  Stranded PyRanges object has 5 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

  >>> ag.sort('col1')
  +--------------+-----------+-----------+--------------+-----------+
  | Chromosome   |     Start |       End | Strand       |      col1 |
  | (category)   |   (int32) |   (int32) | (category)   |   (int64) |
  |--------------+-----------+-----------+--------------+-----------|
  | chrX         |        55 |        88 | +            |         1 |
  | chrX         |        35 |        45 | +            |         2 |
  | chrX         |        75 |        85 | +            |         2 |
  | chrX         |        65 |        75 | +            |         4 |
  | chrX         |        20 |        30 | +            |         4 |
  +--------------+-----------+-----------+--------------+-----------+
  Stranded PyRanges object has 5 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

  >>> ag.sort(['col1', 'End'])
  +--------------+-----------+-----------+--------------+-----------+
  | Chromosome   |     Start |       End | Strand       |      col1 |
  | (category)   |   (int32) |   (int32) | (category)   |   (int64) |
  |--------------+-----------+-----------+--------------+-----------|
  | chrX         |        55 |        88 | +            |         1 |
  | chrX         |        35 |        45 | +            |         2 |
  | chrX         |        75 |        85 | +            |         2 |
  | chrX         |        20 |        30 | +            |         4 |
  | chrX         |        65 |        75 | +            |         4 |
  +--------------+-----------+-----------+--------------+-----------+
  Stranded PyRanges object has 5 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

[add note: index are not allowed. Stil, you can use sort to get rows in a certain order]
Operations on coordinates
[change columns as series: p.Start+=1000 ...]
[... however there are more convenient methods: subsequence, spliced_sequence, extend]
[after extend, show genome_bounds]

Operations on metadata columns:
[insert new columns: 1. p.Col1=... or 2. assign method. 3. Assign with multiple ones at once]

Operations on multiple pyranges
[concatenation: use pandas and turn to pyranges]

A common operation on (multiple) pyranges regard overlaps. These are shown in the next page


Overlapping and matching PyRanges
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[present different methods for different aims that all have to do with overlap: merge, cluster, subtract, join, count_overlaps ... . Start with a table summarizing differences: input, output]. 
[add note: pandas merge: different!]

Summary statistics
~~~~~~~~~~~~~~~~~~

[Create count-matrix from multiple PyRanges]
[all stats methods presented briefly]

Computing with PyRanges
~~~~~~~~~~~~~~~~~~~~~~~

[ready made methods should cover most things]
[possibility to chain things to save memory]
[outline strategies for custom methods: apply and similar methods]
[Also cite the simple but not optimal: convert to dataframes / or iterate through groups of same-chrom dataframes]
[multiple cores]

Working at the transcript level
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[spliced_subsequence, subsequence, get_transcript_sequence, 
extend (to be developed with group_by),
boundaries ,
cumsum groupby as example

]


Fetching external gene tracks 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[if pyranges_db is a thing, describe its uses here]


RLEs: run length encodings
~~~~~~~~~~~~~~~~~~~~~~~~~~

[outline as advanced usage. Put everything related to RLEs in a single chapter; keep as last even if you add further chapters]



