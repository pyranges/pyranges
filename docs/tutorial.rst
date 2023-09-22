Tutorial
========

PyRanges objects represent sequence intervals ("ranges") such as genomic regions.
Examples include gene annotations, mapped reads, protein domains, and results from
Blast or analogous software. PyRanges provides convenient methods to load data in
GFF, GTF, BED, BAM format. In this tutorial, we refer to PyRanges objects for their
typical use case, i.e. representing genomic intervals, but the same concepts and methods
may be applied to RNA or protein sequences.

Every interval in a PyRanges object is defined, at minimum, by its chromosome and start
and end coordinates. Optionally, the strand (+ or -) can also be present; if so, the
PyRanges object is "Stranded". The ranges can additionally have an arbitrary number
of meta-data fields, i.e. columns associated with them.

Note that PyRanges follows the standard python convention:

* coordinates are **0-based**
* in each interval, the **start is included** and the **end is excluded**.

Some genomic coordinate formats (GFF, GTF) use a 1-based, start-and-end-included format.
PyRanges takes care of converting between these conventions when loading and writing files in these formats.
The data in PyRanges objects are stored as `Pandas <https://pandas.pydata.org/>`_ dataframes.
This means that the vast Python ecosystem
for high-performance scientific computing is available to manipulate the data in PyRanges objects.
While not strictly necessary, having familiarity with pandas greatly facilitates the use of PyRanges.
In this tutorial and in how-to pages, we will sometimes mix PyRanges and pandas functionalities.


.. contents:: Contents of Tutorial
   :depth: 3


Getting started
~~~~~~~~~~~~~~~

For this tutorial, we will use real data consisting of the genome sequence (fasta format) and annotation (GFF3 format)
of a worm with a very small but complex genome (*Dimorphilus gyrociliatus*, assembly ID GCA_904063045.1).
These data were modified for this tutorial only to shorten IDs.

Download and unpack tutorial data with:

.. code:: bash

	curl -O https://mariottigenomicslab.bio.ub.edu/pyranges_data/pyranges_tutorial_data.tar.gz
	tar zxf pyranges_tutorial_data.tar.gz
	
We recommend using `ipython <https://ipython.readthedocs.io/>`_ or `Jupyter <https://jupyter.org/>`_ for this tutorial.
Besides pyranges and some of its dependencies, we will use the optional module **pyfaidx**. Make sure to install it with:

.. code-block:: shell

      pip install pyfaidx
      
Or with:

.. code-block:: shell

      conda install -c bioconda pyfaidx
      
Loading and accessing pyranges objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let's import libraries and load the annotation into a PyRanges object called ``ann``:

  >>> import pyranges as pr, pandas as pd
  >>> ann = pr.read_gff3('Dgyro_annotation.gff')
  >>> ann
  +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+--------------------+---------------+------------+-------------+-------+
  | Chromosome   | Source     | Feature      | Start     | End       | Score      | Strand       | Frame      | ID                 | Dbxref        | gbkey      | mol_type    | +13   |
  | (category)   | (object)   | (category)   | (int64)   | (int64)   | (object)   | (category)   | (object)   | (object)           | (object)      | (object)   | (object)    | ...   |
  |--------------+------------+--------------+-----------+-----------+------------+--------------+------------+--------------------+---------------+------------+-------------+-------|
  | 0001.1       | EMBL       | region       | 0         | 3609319   | .          | +            | .          | 0001.1:1..3609319  | taxon:2664684 | Src        | genomic DNA | ...   |
  | 0001.1       | EMBL       | gene         | 7326      | 7606      | .          | +            | .          | gene-DGYR_LOCUS1   | nan           | Gene       | nan         | ...   |
  | 0001.1       | EMBL       | mRNA         | 7326      | 7606      | .          | +            | .          | rna-DGYR_LOCUS1    | nan           | mRNA       | nan         | ...   |
  | 0001.1       | EMBL       | exon         | 7326      | 7606      | .          | +            | .          | exon-DGYR_LOCUS1-1 | nan           | mRNA       | nan         | ...   |
  | ...          | ...        | ...          | ...       | ...       | ...        | ...          | ...        | ...                | ...           | ...        | ...         | ...   |
  | 0346.1       | EMBL       | region       | 0         | 1411      | .          | +            | .          | 0346.1:1..1411     | taxon:2664684 | Src        | genomic DNA | ...   |
  | 0347.1       | EMBL       | region       | 0         | 1094      | .          | +            | .          | 0347.1:1..1094     | taxon:2664684 | Src        | genomic DNA | ...   |
  | 0348.1       | EMBL       | region       | 0         | 667       | .          | +            | .          | 0348.1:1..667      | taxon:2664684 | Src        | genomic DNA | ...   |
  | 0349.1       | EMBL       | region       | 0         | 641       | .          | +            | .          | 0349.1:1..641      | taxon:2664684 | Src        | genomic DNA | ...   |
  +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+--------------------+---------------+------------+-------------+-------+
  Stranded PyRanges object has 241,137 rows and 25 columns from 349 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  13 hidden columns: note, Name, gene_biotype, locus_tag, Parent, Note, standard_name, product, protein_id, pseudo, partial, start_range, end_range


The ``ann`` object has lots of columns, most of which are not displayed because of space. Here's the full list:

  >>> ann.columns
  Index(['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand',
         'Frame', 'ID', 'Dbxref', 'gbkey', 'mol_type', 'note', 'Name',
         'gene_biotype', 'locus_tag', 'Parent', 'Note', 'standard_name',
         'product', 'protein_id', 'pseudo', 'partial', 'start_range',
         'end_range'],
        dtype='object')
      
      
Let's select only certain columns:

  >>> ann = ann[ ['Feature', 'Parent', 'ID'] ]
  >>> ann
  +--------------+--------------+-----------+-----------+--------------+------------------+--------------------+
  | Chromosome   | Feature      | Start     | End       | Strand       | Parent           | ID                 |
  | (category)   | (category)   | (int64)   | (int64)   | (category)   | (object)         | (object)           |
  |--------------+--------------+-----------+-----------+--------------+------------------+--------------------|
  | 0001.1       | region       | 0         | 3609319   | +            | nan              | 0001.1:1..3609319  |
  | 0001.1       | gene         | 7326      | 7606      | +            | nan              | gene-DGYR_LOCUS1   |
  | 0001.1       | mRNA         | 7326      | 7606      | +            | gene-DGYR_LOCUS1 | rna-DGYR_LOCUS1    |
  | 0001.1       | exon         | 7326      | 7606      | +            | rna-DGYR_LOCUS1  | exon-DGYR_LOCUS1-1 |
  | ...          | ...          | ...       | ...       | ...          | ...              | ...                |
  | 0346.1       | region       | 0         | 1411      | +            | nan              | 0346.1:1..1411     |
  | 0347.1       | region       | 0         | 1094      | +            | nan              | 0347.1:1..1094     |
  | 0348.1       | region       | 0         | 667       | +            | nan              | 0348.1:1..667      |
  | 0349.1       | region       | 0         | 641       | +            | nan              | 0349.1:1..641      |
  +--------------+--------------+-----------+-----------+--------------+------------------+--------------------+
  Stranded PyRanges object has 241,137 rows and 7 columns from 349 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

As seen above, column selection syntax is analogous to pandas.
However, a difference is that PyRanges retained the essential columns **Chromosome, Start, End, Strand** even though we did not select them.

The Chromosome column can take any value among the sequence names in the genome assembly.
Only for the best quality assemblies it corresponds to actual chromosomes, and in other cases it is contigs or scaffolds;
for simplicity, here we refer to it as chromosomes. In a fasta file, the sequence name is the first word of a header line
(i.e. those starting with ">"). We can have a peek in the assembly in bash:

.. code:: bash

	grep ">" Dgyro_genome.fa | head
	>0001.1 Dimorphilus gyrociliatus genome assembly, contig: scaffold001, whole genome shotgun sequence
	>0002.1 Dimorphilus gyrociliatus genome assembly, contig: scaffold002, whole genome shotgun sequence
	>0003.1 Dimorphilus gyrociliatus genome assembly, contig: scaffold003, whole genome shotgun sequence
	>0004.1 Dimorphilus gyrociliatus genome assembly, contig: scaffold004, whole genome shotgun sequence
	>0005.1 Dimorphilus gyrociliatus genome assembly, contig: scaffold005, whole genome shotgun sequence
	>0006.1 Dimorphilus gyrociliatus genome assembly, contig: scaffold006, whole genome shotgun sequence
	>0007.1 Dimorphilus gyrociliatus genome assembly, contig: scaffold007, whole genome shotgun sequence
	>0008.1 Dimorphilus gyrociliatus genome assembly, contig: scaffold008, whole genome shotgun sequence
	>0009.1 Dimorphilus gyrociliatus genome assembly, contig: scaffold009, whole genome shotgun sequence
	>0010.1 Dimorphilus gyrociliatus genome assembly, contig: scaffold010, whole genome shotgun sequence

	
Genomic annotations often contain information for diverse entities, such as genes, mRNAs, exons, CDS, etc.
In GFF files, the entity type is encoded in the Feature column. In pyranges, you use the dot notation to
fetch an individual column, which is technically a pandas Series:

  >>> ann.Feature
  0    region
  1      gene
  2      mRNA
  3      exon
  4       CDS
        ...  
  0    region
  0    region
  0    region
  0    region
  0    region
  Name: Feature, Length: 241137, dtype: category
  Categories (6, object): ['CDS', 'exon', 'gene', 'mRNA', 'pseudogene', 'region']
  

Let's focus on a subset of the annotation: CDS intervals, corresponding to coding sequences.
We filter rows and create a new PyRanges object called ``cds``:

  >>> selector = (ann.Feature == 'CDS')
  >>> cds = ann [selector]
  

We used another syntax familiar to pandas users. The object ``selector`` is a Series of boolean values, so it can be used to index PyRanges.

Let's further reduce the width of the cds object. We showcase an alternative method for column selection: the method `drop` lets us choose which columns to discard.

  >>> cds = cds.drop( ['Feature', 'Parent'] )
  >>> cds
  +--------------+-----------+-----------+--------------+------------------+
  | Chromosome   | Start     | End       | Strand       | ID               |
  | (category)   | (int64)   | (int64)   | (category)   | (object)         |
  |--------------+-----------+-----------+--------------+------------------|
  | 0001.1       | 7327      | 7606      | +            | cds-CAD5110615.1 |
  | 0001.1       | 46377     | 47603     | +            | cds-CAD5110625.1 |
  | 0001.1       | 48406     | 49448     | +            | cds-CAD5110625.1 |
  | 0001.1       | 46377     | 47603     | +            | cds-CAD5110626.1 |
  | ...          | ...       | ...       | ...          | ...              |
  | 0117.1       | 26816     | 27881     | -            | cds-CAD5126988.1 |
  | 0117.1       | 36183     | 38697     | -            | cds-CAD5126989.1 |
  | 0117.1       | 39309     | 39450     | -            | cds-CAD5126990.1 |
  | 0117.1       | 38911     | 39256     | -            | cds-CAD5126990.1 |
  +--------------+-----------+-----------+--------------+------------------+
  Stranded PyRanges object has 100,040 rows and 5 columns from 117 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

  
The interface shown so far is analogous to pandas.
Additionally, pyranges offers a non-pandas syntax for selecting intervals in a genomic region of interest (i.e. region retrieval).
The code below will show only intervals completely included in the specified position range in the requested chromosome:

  >>> cds['0002.1', 145000:150000]
  +--------------+-----------+-----------+--------------+------------------+
  |   Chromosome |     Start |       End | Strand       | ID               |
  |   (category) |   (int64) |   (int64) | (category)   | (object)         |
  |--------------+-----------+-----------+--------------+------------------|
  |          2.1 |    146123 |    146246 | +            | cds-CAD5111630.1 |
  |          2.1 |    146306 |    147077 | +            | cds-CAD5111630.1 |
  |          2.1 |    147830 |    147976 | +            | cds-CAD5111631.1 |
  |          2.1 |    148191 |    148297 | +            | cds-CAD5111631.1 |
  |          2.1 |    148360 |    148489 | +            | cds-CAD5111631.1 |
  |          2.1 |    145002 |    145116 | -            | cds-CAD5111629.1 |
  |          2.1 |    145176 |    145262 | -            | cds-CAD5111629.1 |
  |          2.1 |    145320 |    145435 | -            | cds-CAD5111629.1 |
  +--------------+-----------+-----------+--------------+------------------+
  Stranded PyRanges object has 8 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

  
The syntax for region retrieval may consists of:

* chromosome
* chromosome, position slice (as the example above)
* chromosome, strand, position slice

So, for example, this is also valid:

  >>> cds['0002.1', "-", 145000:150000]
  +--------------+-----------+-----------+--------------+------------------+
  |   Chromosome |     Start |       End | Strand       | ID               |
  |   (category) |   (int64) |   (int64) | (category)   | (object)         |
  |--------------+-----------+-----------+--------------+------------------|
  |          2.1 |    145002 |    145116 | -            | cds-CAD5111629.1 |
  |          2.1 |    145176 |    145262 | -            | cds-CAD5111629.1 |
  |          2.1 |    145320 |    145435 | -            | cds-CAD5111629.1 |
  +--------------+-----------+-----------+--------------+------------------+
  Stranded PyRanges object has 3 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

  
It is important to differentiate between **Stranded and Unstranded** PyRanges objects.
When a Strand column is present and all its values are "+" or "-", the object is Stranded.
When there are invalid values (e.g. ".") or the Strand column is absent, it is Unstranded.
You can check whether an interval is Stranded with:

  >>> cds.stranded
  True
  

Certain pyranges methods require a Stranded input.
While the annotation used in this tutorial is naturally Stranded, others may not be.
If necessary, you may use method ``make_stranded`` to transform all invalid Strand values to "+" or remove them.

Working with groups of exons
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Multi-exonic genes are represented with multiple rows in PyRanges. In this tutorial, the ``ID`` column links the
intervals belonging to the same CDS: these rows have the same ID value.
While this concept applies to all annotations, files from different sources may use different column names for this purpose (e.g. transcript_id).
Note that here we focus on CDS regions. These may encompass multiple exons, but they do not span the whole mRNA: the 5'UTRs and 3'UTRs are not included.

Next, we will examine the first and last codon of annotated CDSs. We will obtain their genomic coordinate, then fetch their sequence. 

Method ``spliced_subsequence`` allows to obtain a subregion of groups of intervals. The code below derives the first codon of each CDS group (grouping is defined by their ID):

  >>> first=cds.spliced_subsequence(start=0, end=3, by='ID')
  >>> first
  +--------------+-----------+-----------+--------------+------------------+
  | Chromosome   | Start     | End       | Strand       | ID               |
  | (category)   | (int64)   | (int64)   | (category)   | (object)         |
  |--------------+-----------+-----------+--------------+------------------|
  | 0001.1       | 7327      | 7330      | +            | cds-CAD5110615.1 |
  | 0001.1       | 46377     | 46380     | +            | cds-CAD5110625.1 |
  | 0001.1       | 46377     | 46380     | +            | cds-CAD5110626.1 |
  | 0001.1       | 60099     | 60102     | +            | cds-CAD5110627.1 |
  | ...          | ...       | ...       | ...          | ...              |
  | 0117.1       | 38694     | 38697     | -            | cds-CAD5126989.1 |
  | 0117.1       | 27878     | 27881     | -            | cds-CAD5126988.1 |
  | 0117.1       | 21258     | 21261     | -            | cds-CAD5126985.1 |
  | 0117.1       | 14274     | 14277     | -            | cds-CAD5126984.1 |
  +--------------+-----------+-----------+--------------+------------------+
  Stranded PyRanges object has 16,391 rows and 5 columns from 117 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.


Let's **fetch the sequence** for each of these intervals from our genome fasta file.
The function ``get_sequence`` returns one sequence per interval, which we assign to a new column of our pyranges object:

  >>> first.Sequence = pr.get_sequence(first, 'Dgyro_genome.fa')
  >>> first
  +--------------+-----------+-----------+--------------+------------------+------------+
  | Chromosome   | Start     | End       | Strand       | ID               | Sequence   |
  | (category)   | (int64)   | (int64)   | (category)   | (object)         | (object)   |
  |--------------+-----------+-----------+--------------+------------------+------------|
  | 0001.1       | 7327      | 7330      | +            | cds-CAD5110615.1 | ATG        |
  | 0001.1       | 46377     | 46380     | +            | cds-CAD5110625.1 | ATG        |
  | 0001.1       | 46377     | 46380     | +            | cds-CAD5110626.1 | ATG        |
  | 0001.1       | 60099     | 60102     | +            | cds-CAD5110627.1 | ATG        |
  | ...          | ...       | ...       | ...          | ...              | ...        |
  | 0117.1       | 38694     | 38697     | -            | cds-CAD5126989.1 | ATT        |
  | 0117.1       | 27878     | 27881     | -            | cds-CAD5126988.1 | ATG        |
  | 0117.1       | 21258     | 21261     | -            | cds-CAD5126985.1 | ATG        |
  | 0117.1       | 14274     | 14277     | -            | cds-CAD5126984.1 | ATG        |
  +--------------+-----------+-----------+--------------+------------------+------------+
  Stranded PyRanges object has 16,391 rows and 6 columns from 117 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  
  
The ``Sequence`` column is a pandas Series containing strings. We see that the starting codon is ATG in most cases, as expected.
When we check the length of the sequences, we notice that some are not 3-letter long:

  >>> (first.Sequence.str.len() == 3 ).all()
  False
  
  
Let's look at those sequences, using a row selector as before:

  >>> first [ first.Sequence.str.len() != 3 ]
  +--------------+-----------+-----------+--------------+------------------+------------+
  | Chromosome   | Start     | End       | Strand       | ID               | Sequence   |
  | (category)   | (int64)   | (int64)   | (category)   | (object)         | (object)   |
  |--------------+-----------+-----------+--------------+------------------+------------|
  | 0001.1       | 667699    | 667700    | -            | cds-CAD5110773.1 | a          |
  | 0001.1       | 667641    | 667643    | -            | cds-CAD5110773.1 | TG         |
  | 0002.1       | 2632107   | 2632109   | -            | cds-CAD5112293.1 | AT         |
  | 0002.1       | 2631440   | 2631441   | -            | cds-CAD5112293.1 | G          |
  | ...          | ...       | ...       | ...          | ...              | ...        |
  | 0024.1       | 1091339   | 1091341   | -            | cds-CAD5125104.1 | AT         |
  | 0024.1       | 1091163   | 1091164   | -            | cds-CAD5125104.1 | G          |
  | 0025.1       | 39753     | 39755     | -            | cds-CAD5125115.1 | at         |
  | 0025.1       | 39692     | 39693     | -            | cds-CAD5125115.1 | g          |
  +--------------+-----------+-----------+--------------+------------------+------------+
  Stranded PyRanges object has 26 rows and 6 columns from 11 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  

In some cases the starting codon is split between two exons.
This is uncommon, but we are looking at all protein coding genes of a species, so it is expected at least in a few cases.
How do we get the full codon sequence?

Instead of ``get_sequence``, let's use ``get_transcript_sequence``, which returns the concatenated sequence of a group of intervals,
i.e. joining exons together. The sequence is given 5' to 3'.

  >>> seq_first = pr.get_transcript_sequence(
  ...      first, 
  ...      group_by='ID',
  ...      path='Dgyro_genome.fa'
  ...     )
  >>> seq_first
                         ID Sequence
  0        cds-CAD5110614.1      ATG
  1        cds-CAD5110615.1      ATG
  2        cds-CAD5110616.1      atg
  3        cds-CAD5110617.1      atg
  4        cds-CAD5110618.1      ATG
  ...                   ...      ...
  16373  cds-DGYR_LOCUS9540      atg
  16374  cds-DGYR_LOCUS9596      ATG
  16375  cds-DGYR_LOCUS9732      ATG
  16376   cds-DGYR_LOCUS980      ATG
  16377  cds-DGYR_LOCUS9980      ATG
  <BLANKLINE>
  [16378 rows x 2 columns]

  

``seq_first`` is not a PyRanges object, but a pandas DataFrame. It has a column for the group (ID) and one for the sequence.
Here we confirm the sequence length is always 3:


  >>> (seq_first.Sequence.str.len()==3).all()
  True


Finally, let's quantify how many start codons are ATG, using a bit of pandas magic.
First, we make sure the whole sequence is in uppercase characters.
Then, we make a boolean Series ``is_atg`` which has True corresponding to the ATG codons,
then we sum its values to count the instances of True, creating variable ``n_atg``.
We also store the IDs of the CDSs with ATG codons in the variable ``is_atg_ids``. Finally, we print a summary:

  >>> seq_first.Sequence = seq_first.Sequence.str.upper()
  >>> is_atg = (seq_first.Sequence == 'ATG')
  >>> is_atg_ids = seq_first[is_atg].ID
  >>> n_atg = is_atg.sum()
  >>> print(f'There are {n_atg} ATG start codons out of '
  ... f'{len(seq_first)} CDSs => {n_atg/len(seq_first):.2%}')
  There are 16339 ATG start codons out of 16378 CDSs => 99.76%
      

Now, we want to perform an analogous analysis with stop codons. 
First, we get the a pyranges object of the last codon of each CDS.
Conveniently, the method ``spliced_subsequence`` accepts negative arguments to count from the 3',
so we can obtain the last three nucleotides of CDSs with:

  >>> last = cds.spliced_subsequence(start=-3, by='ID')
  

By not providing an ``end`` argument, we requested intervals that reach the very end of each CDS group.
Let's get their sequence as before, then use pandas function ``value_counts`` to count them:

  >>> seq_last = pr.get_transcript_sequence(last, 'ID',
  ...      'Dgyro_genome.fa')
  >>> seq_last.Sequence = seq_last.Sequence.str.upper()
  >>> seq_last.Sequence.value_counts()
  TAA    8986
  TGA    3859
  TAG    3488
  AAA       4
  CTT       3
  AAG       3
  TTA       3
  AGG       2
  GAG       2
  ATA       2
  GTT       2
  CGA       2
  ACA       2
  TAT       2
  CCC       1
  AGT       1
  TCA       1
  TTC       1
  TTT       1
  AGC       1
  CAA       1
  TTG       1
  AAT       1
  GTG       1
  CCA       1
  AAC       1
  TCT       1
  GGC       1
  GAA       1
  CAG       1
  ATG       1
  CTG       1
  Name: Sequence, dtype: int64
  

The canonical stop codons account for the great majority of cases, but there are some other values.
This may warrant a further look into these CDSs. In this tutorial, we'll simply exclude them from our next steps.

Let's gather the IDs of CDSs with a canonical stop, to be used further on:

  >>> is_stop = seq_last.Sequence.isin( {'TAG', 'TAA', 'TGA'} )
  >>> is_stop_ids = seq_last[is_stop].ID
  

Writing coordinates and sequences to the disk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We want to get a "clean" annotation consisting only of canonical CDSs, with an ATG starting codon
and a TAA/TAG/TGA stop codon. First, we put together the IDs of CDSs with these characteristics:

  >>> clean_ids = set(is_atg_ids).intersection(set(is_stop_ids))
  
  
Then we subset the ann pyranges object:

  >>> clean_ann = ann[ann.ID.isin(clean_ids)]
  >>> clean_ann
  +--------------+--------------+-----------+-----------+--------------+---------------------+------------------+
  | Chromosome   | Feature      | Start     | End       | Strand       | Parent              | ID               |
  | (category)   | (category)   | (int64)   | (int64)   | (category)   | (object)            | (object)         |
  |--------------+--------------+-----------+-----------+--------------+---------------------+------------------|
  | 0001.1       | CDS          | 7327      | 7606      | +            | rna-DGYR_LOCUS1     | cds-CAD5110615.1 |
  | 0001.1       | CDS          | 46377     | 47603     | +            | rna-DGYR_LOCUS4     | cds-CAD5110625.1 |
  | 0001.1       | CDS          | 48406     | 49448     | +            | rna-DGYR_LOCUS4     | cds-CAD5110625.1 |
  | 0001.1       | CDS          | 46377     | 47603     | +            | rna-DGYR_LOCUS4-2   | cds-CAD5110626.1 |
  | ...          | ...          | ...       | ...       | ...          | ...                 | ...              |
  | 0117.1       | CDS          | 20172     | 21261     | -            | rna-DGYR_LOCUS14199 | cds-CAD5126985.1 |
  | 0117.1       | CDS          | 26816     | 27881     | -            | rna-DGYR_LOCUS14201 | cds-CAD5126988.1 |
  | 0117.1       | CDS          | 39309     | 39450     | -            | rna-DGYR_LOCUS14203 | cds-CAD5126990.1 |
  | 0117.1       | CDS          | 38911     | 39256     | -            | rna-DGYR_LOCUS14203 | cds-CAD5126990.1 |
  +--------------+--------------+-----------+-----------+--------------+---------------------+------------------+
  Stranded PyRanges object has 99,155 rows and 7 columns from 117 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  

We can now write this pyranges object to a file, for example in GTF format:

  >>> clean_ann.to_gtf('Dgyro_annotation.canonical_CDS.gtf')
  

Let's get the sequence for the canonical CDSs and write it to a tabular file. 

  >>> clean_ann_seq = pr.get_transcript_sequence(clean_ann, 'ID',
  ...           'Dgyro_genome.fa')
  >>> clean_ann_seq.to_csv('Dgyro_canonical_CDS.seq.tsv', 
  ...                sep='\t', index=False)
                     


Note that ``clean_ann_seq`` is a pandas DataFrame. To write sequences in fasta format we use: 

  >>> with open('Dgyro_canonical_CDS.fa', 'w') as fw: # doctest: +SKIP
  ...   for xin, xid, xseq in clean_ann_seq.itertuples():
  ...     fw.write(f'>{xid}\n{xseq}\n')
  
  
Extending genomic intervals
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now we want to obtain (a practical approximation of) promoter sequences, here defined as the 300bp region before the start codon.
Before we begin, let's peek into our object ``cds``:

  >>> cds.head()
  +--------------+-----------+-----------+--------------+------------------+
  |   Chromosome |     Start |       End | Strand       | ID               |
  |   (category) |   (int64) |   (int64) | (category)   | (object)         |
  |--------------+-----------+-----------+--------------+------------------|
  |          1.1 |      7327 |      7606 | +            | cds-CAD5110615.1 |
  |          1.1 |     46377 |     47603 | +            | cds-CAD5110625.1 |
  |          1.1 |     48406 |     49448 | +            | cds-CAD5110625.1 |
  |          1.1 |     46377 |     47603 | +            | cds-CAD5110626.1 |
  |          1.1 |     48406 |     48736 | +            | cds-CAD5110626.1 |
  |          1.1 |     48839 |     48912 | +            | cds-CAD5110626.1 |
  |          1.1 |     60099 |     60409 | +            | cds-CAD5110627.1 |
  |          1.1 |     60476 |     60515 | +            | cds-CAD5110627.1 |
  +--------------+-----------+-----------+--------------+------------------+
  Stranded PyRanges object has 8 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  
  
First, we use the method ``extend`` to obtain intervals which include the CDS and the promoter defined as above:

  >>> g = cds.extend({'5':300}, group_by='ID')
  >>> g.head()
  +--------------+-----------+-----------+--------------+------------------+
  |   Chromosome |     Start |       End | Strand       | ID               |
  |   (category) |   (int64) |   (int64) | (category)   | (object)         |
  |--------------+-----------+-----------+--------------+------------------|
  |          1.1 |      7027 |      7606 | +            | cds-CAD5110615.1 |
  |          1.1 |     46077 |     47603 | +            | cds-CAD5110625.1 |
  |          1.1 |     48406 |     49448 | +            | cds-CAD5110625.1 |
  |          1.1 |     46077 |     47603 | +            | cds-CAD5110626.1 |
  |          1.1 |     48406 |     48736 | +            | cds-CAD5110626.1 |
  |          1.1 |     48839 |     48912 | +            | cds-CAD5110626.1 |
  |          1.1 |     59799 |     60409 | +            | cds-CAD5110627.1 |
  |          1.1 |     60476 |     60515 | +            | cds-CAD5110627.1 |
  +--------------+-----------+-----------+--------------+------------------+
  Stranded PyRanges object has 8 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.


The first argument ensures that the 300bp extension is applied only at the 5' (left side for + strand intervals, right side for - strand intervals).
Through the ``group_by`` argument, we request one extension per CDS, instead of extending every interval.
In the object we obtained, the promoter corresponds to the first 300 bp of every interval group.
We can use method ``spliced_subsequence`` again to get it:

  >>> prom = g.spliced_subsequence(0, 300, 'ID')
  >>> prom.head()
  +--------------+-----------+-----------+--------------+------------------+
  |   Chromosome |     Start |       End | Strand       | ID               |
  |   (category) |   (int64) |   (int64) | (category)   | (object)         |
  |--------------+-----------+-----------+--------------+------------------|
  |          1.1 |      7027 |      7327 | +            | cds-CAD5110615.1 |
  |          1.1 |     46077 |     46377 | +            | cds-CAD5110625.1 |
  |          1.1 |     46077 |     46377 | +            | cds-CAD5110626.1 |
  |          1.1 |     59799 |     60099 | +            | cds-CAD5110627.1 |
  |          1.1 |     59799 |     60099 | +            | cds-CAD5110628.1 |
  |          1.1 |     72099 |     72399 | +            | cds-CAD5110629.1 |
  |          1.1 |     75736 |     76036 | +            | cds-CAD5110630.1 |
  |          1.1 |     79997 |     80297 | +            | cds-CAD5110631.1 |
  +--------------+-----------+-----------+--------------+------------------+
  Stranded PyRanges object has 8 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  

Because we extended intervals, some may have gone out-of-bounds on the left or on the right side:
they may have a Start smaller than 0, or an End greater than the length of its chromosome, respectively.
Indeed, we see there are cases of the first type:

  >>> prom[prom.Start<0]
  +--------------+-----------+-----------+--------------+---------------------+
  |   Chromosome |     Start |       End | Strand       | ID                  |
  |   (category) |   (int64) |   (int64) | (category)   | (object)            |
  |--------------+-----------+-----------+--------------+---------------------|
  |         12.1 |      -129 |       171 | +            | cds-DGYR_LOCUS8189  |
  |         18.1 |       -12 |       288 | +            | cds-DGYR_LOCUS10365 |
  |         48.1 |      -118 |       182 | +            | cds-CAD5126431.1    |
  |         78.1 |      -299 |         1 | +            | cds-CAD5126743.1    |
  +--------------+-----------+-----------+--------------+---------------------+
  Stranded PyRanges object has 4 rows and 5 columns from 4 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.


The function ``genome_bounds`` in submodule ``genomicfeatures`` is designed to correct this.
We may use it to remove out-of-bounds intervals, or to retain only their in-bound portions. We go for the second option, with ``clip=True``:

  >>> import pyfaidx
  >>> pyf=pyfaidx.Fasta('Dgyro_genome.fa')
  >>> cor_prom = pr.genomicfeatures.genome_bounds(prom,
  ...               chromsizes=pyf,
  ...               clip=True)
                    
                    

To detect cases of out-of-bounds on the right side, function ``genome_bounds`` needs to know chromosome sizes.
Various input types are accepted for the ``chromsizes`` argument; we used a ``pyfaidx.Fasta`` object, which derives it from a fasta file.

The intervals above (and also the right-side out-of-bounds, though we don't inspect them) have been corrected:


  >>> outofbounds_left=prom[prom.Start<0].ID
  >>> cor_prom[cor_prom.ID.isin(outofbounds_left)]
  +--------------+-----------+-----------+--------------+---------------------+
  |   Chromosome |     Start |       End | Strand       | ID                  |
  |   (category) |   (int64) |   (int64) | (category)   | (object)            |
  |--------------+-----------+-----------+--------------+---------------------|
  |         12.1 |         0 |       171 | +            | cds-DGYR_LOCUS8189  |
  |         18.1 |         0 |       288 | +            | cds-DGYR_LOCUS10365 |
  |         48.1 |         0 |       182 | +            | cds-CAD5126431.1    |
  |         78.1 |         0 |         1 | +            | cds-CAD5126743.1    |
  +--------------+-----------+-----------+--------------+---------------------+
  Stranded PyRanges object has 4 rows and 5 columns from 4 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.


Detecting overlaps among intervals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let's see if any of the promoter regions overlap other CDSs. Pyranges offers many efficient methods to detect overlaps, such as ``overlap``:

  >>> cor_prom.overlap(cds, strandedness=False)
  +--------------+-----------+-----------+--------------+------------------+
  | Chromosome   | Start     | End       | Strand       | ID               |
  | (category)   | (int64)   | (int64)   | (category)   | (object)         |
  |--------------+-----------+-----------+--------------+------------------|
  | 0001.1       | 320232    | 320532    | +            | cds-CAD5110677.1 |
  | 0001.1       | 371611    | 371911    | +            | cds-CAD5110692.1 |
  | 0001.1       | 434425    | 434725    | +            | cds-CAD5110703.1 |
  | 0001.1       | 445177    | 445477    | +            | cds-CAD5110709.1 |
  | ...          | ...       | ...       | ...          | ...              |
  | 0111.1       | 42019     | 42319     | -            | cds-CAD5126964.1 |
  | 0114.1       | 4745      | 5045      | -            | cds-CAD5126973.1 |
  | 0117.1       | 12844     | 13144     | +            | cds-CAD5126983.1 |
  | 0117.1       | 38697     | 38997     | -            | cds-CAD5126989.1 |
  +--------------+-----------+-----------+--------------+------------------+
  Stranded PyRanges object has 2,374 rows and 5 columns from 61 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  
  
By default, this method reports intervals in the self pyranges object (i.e., ``cor_prom``) that have at least 1bp of overlap with
the other pyranges (i.e., cds). By invoking ``strandedness=False``, we included overlaps even between intervals on opposite strands.

There are many promoters overlapping CDSs. Let's get the overlapping regions only, using function ``intersect``:

  >>> prom_in_cds = cor_prom.intersect(cds, strandedness=False)
  >>> prom_in_cds
  +--------------+-----------+-----------+--------------+------------------+
  | Chromosome   | Start     | End       | Strand       | ID               |
  | (category)   | (int64)   | (int64)   | (category)   | (object)         |
  |--------------+-----------+-----------+--------------+------------------|
  | 0001.1       | 320529    | 320532    | +            | cds-CAD5110677.1 |
  | 0001.1       | 371611    | 371744    | +            | cds-CAD5110692.1 |
  | 0001.1       | 371611    | 371744    | +            | cds-CAD5110692.1 |
  | 0001.1       | 371870    | 371911    | +            | cds-CAD5110692.1 |
  | ...          | ...       | ...       | ...          | ...              |
  | 0114.1       | 4910      | 5045      | -            | cds-CAD5126973.1 |
  | 0117.1       | 13016     | 13019     | +            | cds-CAD5126983.1 |
  | 0117.1       | 13080     | 13144     | +            | cds-CAD5126983.1 |
  | 0117.1       | 38911     | 38997     | -            | cds-CAD5126989.1 |
  +--------------+-----------+-----------+--------------+------------------+
  Stranded PyRanges object has 4,185 rows and 5 columns from 61 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  
  
``intersect`` returned more rows than ``overlap``. This is because an interval in ``cor_prom`` may overlap multiple intervals in cds,
potentially with different intersection regions (compare the 2nd, 3rd and 4th rows, which are all subregions of the same promoter).
Therefore, ``intersect`` returns one row for each pair of overlapping intervals, while `overlap` always returns a subset of rows from the self pyranges object, unaltered.

We want to remove redundancy in the object above. We use the method ``merge`` to fuse intervals that have some overlap (and the same ID value):

  >>> prom_in_cds = prom_in_cds.merge(by='ID')
  >>> prom_in_cds
  +--------------+-----------+-----------+--------------+------------------+
  | Chromosome   | Start     | End       | Strand       | ID               |
  | (category)   | (int64)   | (int64)   | (category)   | (object)         |
  |--------------+-----------+-----------+--------------+------------------|
  | 0001.1       | 320529    | 320532    | +            | cds-CAD5110677.1 |
  | 0001.1       | 371611    | 371744    | +            | cds-CAD5110692.1 |
  | 0001.1       | 371870    | 371911    | +            | cds-CAD5110692.1 |
  | 0001.1       | 434527    | 434725    | +            | cds-CAD5110703.1 |
  | ...          | ...       | ...       | ...          | ...              |
  | 0114.1       | 4910      | 5045      | -            | cds-CAD5126973.1 |
  | 0117.1       | 13016     | 13019     | +            | cds-CAD5126983.1 |
  | 0117.1       | 13080     | 13144     | +            | cds-CAD5126983.1 |
  | 0117.1       | 38911     | 38997     | -            | cds-CAD5126989.1 |
  +--------------+-----------+-----------+--------------+------------------+
  Stranded PyRanges object has 3,040 rows and 5 columns from 61 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

  
Note that with overlap or intersect, we do not keep track of the coordinates of overlapping intervals in the second object (``cds``),
as we only obtain those in the first object (``cor_prom``). For that task, check methods ``cluster`` and ``join`` (not shown here).

We now want to calculate how long the overlapping regions are. We create a new column named ``Length`` by subtracting ``Start`` from ``End``.
This operation is performed element-wise (the 1st value of End minus the 1st value of Start, the 2nd End minus the 2nd Start, etc), a common paradigm of pandas Series.

  >>> prom_in_cds.Length = prom_in_cds.End - prom_in_cds.Start
  >>> prom_in_cds
  +--------------+-----------+-----------+--------------+------------------+-----------+
  | Chromosome   | Start     | End       | Strand       | ID               | Length    |
  | (category)   | (int64)   | (int64)   | (category)   | (object)         | (int64)   |
  |--------------+-----------+-----------+--------------+------------------+-----------|
  | 0001.1       | 320529    | 320532    | +            | cds-CAD5110677.1 | 3         |
  | 0001.1       | 371611    | 371744    | +            | cds-CAD5110692.1 | 133       |
  | 0001.1       | 371870    | 371911    | +            | cds-CAD5110692.1 | 41        |
  | 0001.1       | 434527    | 434725    | +            | cds-CAD5110703.1 | 198       |
  | ...          | ...       | ...       | ...          | ...              | ...       |
  | 0114.1       | 4910      | 5045      | -            | cds-CAD5126973.1 | 135       |
  | 0117.1       | 13016     | 13019     | +            | cds-CAD5126983.1 | 3         |
  | 0117.1       | 13080     | 13144     | +            | cds-CAD5126983.1 | 64        |
  | 0117.1       | 38911     | 38997     | -            | cds-CAD5126989.1 | 86        |
  +--------------+-----------+-----------+--------------+------------------+-----------+
  Stranded PyRanges object has 3,040 rows and 6 columns from 61 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

  
Pandas vs Pyranges
~~~~~~~~~~~~~~~~~~
It is convenient to think of PyRanges objects as pandas DataFrames decorated with convenient methods for genomic analyses.
As seen above, PyRanges offers an interface analogous to DataFrame for data access and input/output, and it is similar when printed.
Also, the columns of both object types are pandas Series.

Yet, PyRanges is not implemented as a subclass of DataFrame, as we will see shortly, so that it does not offer all its methods.
When in need of a pandas functionality missing in PyRanges, you can easily obtain a DataFrame version of it with property ``df`` (a shortcut for method ``as_df``).
Note that this copies all data: avoid it if you can stick to PyRanges functions.

  >>> prom_in_cds.df
       Chromosome   Start     End Strand                ID  Length
  0        0001.1  320529  320532      +  cds-CAD5110677.1       3
  1        0001.1  371611  371744      +  cds-CAD5110692.1     133
  2        0001.1  371870  371911      +  cds-CAD5110692.1      41
  3        0001.1  434527  434725      +  cds-CAD5110703.1     198
  4        0001.1  445384  445477      +  cds-CAD5110709.1      93
  ...         ...     ...     ...    ...               ...     ...
  3035     0111.1   42092   42201      -  cds-CAD5126964.1     109
  3036     0114.1    4910    5045      -  cds-CAD5126973.1     135
  3037     0117.1   13016   13019      +  cds-CAD5126983.1       3
  3038     0117.1   13080   13144      +  cds-CAD5126983.1      64
  3039     0117.1   38911   38997      -  cds-CAD5126989.1      86
  <BLANKLINE>
  [3040 rows x 6 columns]
  

Let's use the DataFrame for a ``groupby`` operation wherein we get the aggregated length per promoter of regions overlapping a CDS, as pandas Series:

  >>> tot_len = prom_in_cds.df.groupby("ID").Length.sum()
  >>> tot_len.name = 'Tot_length'
  >>> tot_len
  ID
  cds-CAD5110617.1         73
  cds-CAD5110618.1        235
  cds-CAD5110619.1        235
  cds-CAD5110622.1         73
  cds-CAD5110623.1         73
                         ... 
  cds-DGYR_LOCUS5056       14
  cds-DGYR_LOCUS5675       41
  cds-DGYR_LOCUS7571-2    234
  cds-DGYR_LOCUS9062        4
  cds-DGYR_LOCUS980        99
  Name: Tot_length, Length: 2374, dtype: int64
  
  
Let's add this new information (how much of a CDS promoter is overlapping a different CDS) to the ``cds`` object.
Since it is one number per CDS, all intervals with the same ID will have the same ``Tot_length``. This operation corresponds to a database "join",
which is missing from PyRanges functionalities but available as pandas ``merge``:

  >>> z = cds.df.merge(tot_len, on='ID', how='left')
  >>> z.Tot_length.fillna(0, inplace=True, downcast='infer')
  >>> z
         Chromosome  Start    End Strand                ID  Tot_length
  0          0001.1   7327   7606      +  cds-CAD5110615.1           0
  1          0001.1  46377  47603      +  cds-CAD5110625.1           0
  2          0001.1  48406  49448      +  cds-CAD5110625.1           0
  3          0001.1  46377  47603      +  cds-CAD5110626.1           0
  4          0001.1  48406  48736      +  cds-CAD5110626.1           0
  ...           ...    ...    ...    ...               ...         ...
  100035     0117.1  20172  21261      -  cds-CAD5126985.1           0
  100036     0117.1  26816  27881      -  cds-CAD5126988.1           0
  100037     0117.1  36183  38697      -  cds-CAD5126989.1          86
  100038     0117.1  39309  39450      -  cds-CAD5126990.1           0
  100039     0117.1  38911  39256      -  cds-CAD5126990.1           0
  <BLANKLINE>
  [100040 rows x 6 columns]
  
  
Only some CDSs have a promoter overlapping another CDS, so we used how='left' when calling ``merge``.
This retains all rows of ``cds``, introducing NaN values for IDs missing in `tot_len`. On the next line of code, we filled NaN with zeros.

Now let's convert the resulting DataFrame ``z`` back to PyRanges:

  >>> cds = pr.PyRanges(z)
  >>> cds
  +--------------+-----------+-----------+--------------+------------------+--------------+
  | Chromosome   | Start     | End       | Strand       | ID               | Tot_length   |
  | (category)   | (int64)   | (int64)   | (category)   | (object)         | (int64)      |
  |--------------+-----------+-----------+--------------+------------------+--------------|
  | 0001.1       | 7327      | 7606      | +            | cds-CAD5110615.1 | 0            |
  | 0001.1       | 46377     | 47603     | +            | cds-CAD5110625.1 | 0            |
  | 0001.1       | 48406     | 49448     | +            | cds-CAD5110625.1 | 0            |
  | 0001.1       | 46377     | 47603     | +            | cds-CAD5110626.1 | 0            |
  | ...          | ...       | ...       | ...          | ...              | ...          |
  | 0117.1       | 26816     | 27881     | -            | cds-CAD5126988.1 | 0            |
  | 0117.1       | 36183     | 38697     | -            | cds-CAD5126989.1 | 86           |
  | 0117.1       | 39309     | 39450     | -            | cds-CAD5126990.1 | 0            |
  | 0117.1       | 38911     | 39256     | -            | cds-CAD5126990.1 | 0            |
  +--------------+-----------+-----------+--------------+------------------+--------------+
  Stranded PyRanges object has 100,040 rows and 6 columns from 117 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  
  
Now let's dig into the **differences of PyRanges and DataFrame**.
Say that we want to order our intervals by ``Tot_length``. We use PyRanges method ``sort``:


  >>> srt_cds = cds.sort('Tot_length')
  >>> srt_cds
  +--------------+-----------+-----------+--------------+------------------+--------------+
  | Chromosome   | Start     | End       | Strand       | ID               | Tot_length   |
  | (category)   | (int64)   | (int64)   | (category)   | (object)         | (int64)      |
  |--------------+-----------+-----------+--------------+------------------+--------------|
  | 0001.1       | 7327      | 7606      | +            | cds-CAD5110615.1 | 0            |
  | 0001.1       | 2092958   | 2093126   | +            | cds-CAD5111233.1 | 0            |
  | 0001.1       | 2093184   | 2093293   | +            | cds-CAD5111233.1 | 0            |
  | 0001.1       | 2093350   | 2093443   | +            | cds-CAD5111233.1 | 0            |
  | ...          | ...       | ...       | ...          | ...              | ...          |
  | 0117.1       | 26816     | 27881     | -            | cds-CAD5126988.1 | 0            |
  | 0117.1       | 39309     | 39450     | -            | cds-CAD5126990.1 | 0            |
  | 0117.1       | 38911     | 39256     | -            | cds-CAD5126990.1 | 0            |
  | 0117.1       | 36183     | 38697     | -            | cds-CAD5126989.1 | 86           |
  +--------------+-----------+-----------+--------------+------------------+--------------+
  Stranded PyRanges object has 100,040 rows and 6 columns from 117 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  
  
If we sort the analogous DataFrame using pandas ``sort_values``, we see the results does not quite look the same:

  >>> cds.df.sort_values('Tot_length')
        Chromosome    Start      End Strand                ID  Tot_length
  0         0001.1     7327     7606      +  cds-CAD5110615.1           0
  64767     0013.1  1245868  1246381      -  cds-CAD5121003.1           0
  64766     0013.1  1246440  1246707      -  cds-CAD5121003.1           0
  64765     0013.1  1246768  1246871      -  cds-CAD5121003.1           0
  64764     0013.1  1246936  1247072      -  cds-CAD5121003.1           0
  ...          ...      ...      ...    ...               ...         ...
  9380      0002.1   258719   259098      -  cds-CAD5111666.1         300
  64403     0013.1   777119   777332      -  cds-CAD5120886.1         300
  9381      0002.1   258298   258496      -  cds-CAD5111666.1         300
  64404     0013.1   776667   776756      -  cds-CAD5120886.1         300
  48837     0009.1   731837   731964      +  cds-CAD5118662.1         300
  <BLANKLINE>
  [100040 rows x 6 columns]
  
  
Why is that? Under the hood, each PyRanges object is a **collection of DataFrames**: data is spread into separate tables,
one for each chromosome/strand pair; e.g. there is a DataFrame with coordinates for + intervals on chr1, one for - intervals on chr1,
one for + on chr2, one for - on chr2, etc. The user typically does not need to directly access them (but if you do,
you can check the dictionary-type PyRanges attribute ``dfs``).

Intervals on different chromosomes have no order relative to each other, and they are never mixed up in the same table. Indeed, when inspecting a PyRanges object, you see the message::

	For printing, the PyRanges was sorted on Chromosome and Strand.

PyRanges ``sort`` therefore acts on each table independently. Pandas, on the other hand,
has no problems mixing up rows corresponding to different chromosomes, which explains the discrepancy seen above.

This leads to another important difference. In pandas, the index is an essential component of the DataFrame,
providing the row labels, order, and a tool for data access.
**In PyRanges objects, there is no index**. Their internal tables of course have their own indices, but they
are purposely hidden from the user as they are not to be queried or relied upon.

The user should also beware of methods ``merge`` and ``join``, which have different meanings.
As seen above, in Pandas they are slight variations of database "join", while in PyRanges they refer to interval manipulation based on genomic overlap.

Method chaining and custom functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Like Pandas, pyranges support method chaining.
Some useful methods in this sense are: ``pc``, which prints the object and returns it for further chaining;
``subset``, which performs row selection; and `assign`, which adds a column.

We may chain these to obtain at once the CDS subset, add the length of intervals as new column, drop a couple of columns, then print before and after sorting by length:

  >>> ( ann.subset(lambda x:x.Feature=='CDS')
  ... .assign('Length', lambda x:x.End-x.Start)
  ... .drop(['Parent', 'Feature'])
  ... .pc()
  ... .sort('Length')
  ... .print()
  ... )
  +--------------+-----------+-----------+--------------+------------------+-----------+
  | Chromosome   | Start     | End       | Strand       | ID               | Length    |
  | (category)   | (int64)   | (int64)   | (category)   | (object)         | (int64)   |
  |--------------+-----------+-----------+--------------+------------------+-----------|
  | 0001.1       | 7327      | 7606      | +            | cds-CAD5110615.1 | 279       |
  | 0001.1       | 46377     | 47603     | +            | cds-CAD5110625.1 | 1226      |
  | 0001.1       | 48406     | 49448     | +            | cds-CAD5110625.1 | 1042      |
  | 0001.1       | 46377     | 47603     | +            | cds-CAD5110626.1 | 1226      |
  | ...          | ...       | ...       | ...          | ...              | ...       |
  | 0117.1       | 26816     | 27881     | -            | cds-CAD5126988.1 | 1065      |
  | 0117.1       | 36183     | 38697     | -            | cds-CAD5126989.1 | 2514      |
  | 0117.1       | 39309     | 39450     | -            | cds-CAD5126990.1 | 141       |
  | 0117.1       | 38911     | 39256     | -            | cds-CAD5126990.1 | 345       |
  +--------------+-----------+-----------+--------------+------------------+-----------+
  Stranded PyRanges object has 100,040 rows and 6 columns from 117 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  +--------------+-----------+-----------+--------------+------------------+-----------+
  | Chromosome   | Start     | End       | Strand       | ID               | Length    |
  | (category)   | (int64)   | (int64)   | (category)   | (object)         | (int64)   |
  |--------------+-----------+-----------+--------------+------------------+-----------|
  | 0001.1       | 2789047   | 2789050   | +            | cds-CAD5111421.1 | 3         |
  | 0001.1       | 909263    | 909266    | +            | cds-CAD5110841.1 | 3         |
  | 0001.1       | 1721311   | 1721314   | +            | cds-CAD5111112.1 | 3         |
  | 0001.1       | 651283    | 651286    | +            | cds-CAD5110767.1 | 3         |
  | ...          | ...       | ...       | ...          | ...              | ...       |
  | 0117.1       | 38911     | 39256     | -            | cds-CAD5126990.1 | 345       |
  | 0117.1       | 26816     | 27881     | -            | cds-CAD5126988.1 | 1065      |
  | 0117.1       | 20172     | 21261     | -            | cds-CAD5126985.1 | 1089      |
  | 0117.1       | 36183     | 38697     | -            | cds-CAD5126989.1 | 2514      |
  +--------------+-----------+-----------+--------------+------------------+-----------+
  Stranded PyRanges object has 100,040 rows and 6 columns from 117 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.


With ``assign`` and ``subset``, you provide a function which is applied to each DataFrame in the collection.
In both cases, it must return a Series with the same number of rows as the input PyRanges.
For ``subset``, it must be a boolean Series, which is used as row selector. Another similar method called ``apply`` also exists.
Use ``apply`` when your function returns a DataFrame which can be converted to a PyRanges object (i.e. containing Chromosome, Start, End, Strand columns).

In the following code, we select only CDS intervals whose center point is within 50 nucleotides from the middle of a chromosome.
For this, we precompute a table of chromosome sizes from the object `pyf` used earlier, and we use pandas ``merge``, through pyranges ``apply``,
to join this table with each dataframe.


  >>> chromsizes = pd.DataFrame.from_dict(
  ... {'Chromosome':[k for k,v in pyf.items()],
  ...  'chromsize':[len(v) for k,v in pyf.items()]}
  ...    )

  >>> (ann.subset(lambda x:x.Feature=='CDS')
  ... .drop(['Parent', 'Feature', 'ID'])
  ... .apply(lambda x:x.merge(chromsizes, on='Chromosome'))
  ... .assign('midpoint', lambda x:(x.End+x.Start)/2)
  ... .subset(lambda x:abs(x.midpoint - x.chromsize/2)<50)
  ...  )
  +--------------+-----------+-----------+--------------+-------------+-------------+
  | Chromosome   | Start     | End       | Strand       | chromsize   | midpoint    |
  | (object)     | (int64)   | (int64)   | (category)   | (int64)     | (float64)   |
  |--------------+-----------+-----------+--------------+-------------+-------------|
  | 0003.1       | 1616314   | 1616368   | +            | 3232594     | 1616341.0   |
  | 0005.1       | 2273266   | 2273365   | +            | 4546684     | 2273315.5   |
  | 0013.1       | 861154    | 861506    | +            | 1722566     | 861330.0    |
  | 0014.1       | 1120851   | 1120981   | -            | 2241898     | 1120916.0   |
  | ...          | ...       | ...       | ...          | ...         | ...         |
  | 0034.1       | 98757     | 98889     | -            | 197648      | 98823.0     |
  | 0057.1       | 47117     | 47335     | +            | 94438       | 47226.0     |
  | 0065.1       | 81326     | 81695     | +            | 162924      | 81510.5     |
  | 0098.1       | 63943     | 64035     | +            | 127883      | 63989.0     |
  +--------------+-----------+-----------+--------------+-------------+-------------+
  Stranded PyRanges object has 12 rows and 6 columns from 10 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

    
A common operation in pandas is group by then apply, i.e. dividing the table in groups and performing certain operations on each group.
You can do such operations using ``subset``, ``assign``, or ``apply``, depending on what to do with the result.
Remember that the table of each Chromosome/Strand is processed independently.

Let's use this functionality to get the first (5'-most) exon of each CDS group.
We use pyranges ``sort`` with the argument '5', which puts intervals in 5' -> 3' order, then we apply a groupby+apply pandas chain:


  >>> ( ann.subset(lambda x:x.Feature=='CDS').drop(['Parent', 'Feature']).sort('5').apply(lambda x:x.groupby('ID', as_index=False).first()))
  +------------------+--------------+-----------+-----------+--------------+
  | ID               | Chromosome   | Start     | End       | Strand       |
  | (object)         | (category)   | (int64)   | (int64)   | (category)   |
  |------------------+--------------+-----------+-----------+--------------|
  | cds-CAD5110615.1 | 0001.1       | 7327      | 7606      | +            |
  | cds-CAD5110625.1 | 0001.1       | 46377     | 47603     | +            |
  | cds-CAD5110626.1 | 0001.1       | 46377     | 47603     | +            |
  | cds-CAD5110627.1 | 0001.1       | 60099     | 60409     | +            |
  | ...              | ...          | ...       | ...       | ...          |
  | cds-CAD5126985.1 | 0117.1       | 20172     | 21261     | -            |
  | cds-CAD5126988.1 | 0117.1       | 26816     | 27881     | -            |
  | cds-CAD5126989.1 | 0117.1       | 36183     | 38697     | -            |
  | cds-CAD5126990.1 | 0117.1       | 39309     | 39450     | -            |
  +------------------+--------------+-----------+-----------+--------------+
  Stranded PyRanges object has 16,378 rows and 5 columns from 117 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.


This concludes our tutorial. The next pages will delve into pyranges functionalities grouped by topic.

