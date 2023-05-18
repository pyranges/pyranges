"""Module of example data.

See Also
--------

pyranges.random : generate random PyRanges

Examples
--------

>>> pr.data.f1()
+--------------+-----------+-----------+------------+-----------+--------------+
| Chromosome   |     Start |       End | Name       |     Score | Strand       |
| (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
|--------------+-----------+-----------+------------+-----------+--------------|
| chr1         |         3 |         6 | interval1  |         0 | +            |
| chr1         |         8 |         9 | interval3  |         0 | +            |
| chr1         |         5 |         7 | interval2  |         0 | -            |
+--------------+-----------+-----------+------------+-----------+--------------+
Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
For printing, the PyRanges was sorted on Chromosome and Strand.
"""

import pandas as pd
import pkg_resources

import pyranges as pr

__all__ = [
    "f1",
    "f2",
    "chipseq",
    "chipseq_background",
    "aorta",
    "aorta2",
    "ensembl_gtf",
    "gencode_gtf",
    "ucsc_bed",
    "control_bam",
    "cpg",
    "exons",
    "chromsizes",
]


def get_example_path(basename):
    full_path = pkg_resources.resource_filename("pyranges", "example_data/{}".format(basename))

    if full_path.endswith(".bam"):
        # hack to load index too
        pkg_resources.resource_filename("pyranges", "example_data/{}.bai".format(basename))

    return full_path


def aorta():
    """
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # | Chromosome   | Start     | End       | Name       | Score     | Strand       |
    >>> # | (category)   | (int64)   | (int64)   | (object)   | (int64)   | (category)   |
    >>> # |--------------+-----------+-----------+------------+-----------+--------------|
    >>> # | chr1         | 9939      | 10138     | H3K27me3   | 7         | +            |
    >>> # | chr1         | 9953      | 10152     | H3K27me3   | 5         | +            |
    >>> # | chr1         | 10024     | 10223     | H3K27me3   | 1         | +            |
    >>> # | chr1         | 10246     | 10445     | H3K27me3   | 4         | +            |
    >>> # | ...          | ...       | ...       | ...        | ...       | ...          |
    >>> # | chr1         | 9978      | 10177     | H3K27me3   | 7         | -            |
    >>> # | chr1         | 10001     | 10200     | H3K27me3   | 5         | -            |
    >>> # | chr1         | 10127     | 10326     | H3K27me3   | 1         | -            |
    >>> # | chr1         | 10241     | 10440     | H3K27me3   | 6         | -            |
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # Stranded PyRanges object has 11 rows and 6 columns from 1 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("aorta.bed")

    return pr.read_bed(full_path)


def aorta2():
    """
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # | Chromosome   | Start     | End       | Name       | Score     | Strand       |
    >>> # | (category)   | (int64)   | (int64)   | (object)   | (int64)   | (category)   |
    >>> # |--------------+-----------+-----------+------------+-----------+--------------|
    >>> # | chr1         | 10073     | 10272     | Input      | 1         | +            |
    >>> # | chr1         | 10280     | 10479     | Input      | 1         | +            |
    >>> # | chr1         | 16056     | 16255     | Input      | 1         | +            |
    >>> # | chr1         | 16064     | 16263     | Input      | 1         | +            |
    >>> # | ...          | ...       | ...       | ...        | ...       | ...          |
    >>> # | chr1         | 10079     | 10278     | Input      | 1         | -            |
    >>> # | chr1         | 10082     | 10281     | Input      | 1         | -            |
    >>> # | chr1         | 10149     | 10348     | Input      | 1         | -            |
    >>> # | chr1         | 19958     | 20157     | Input      | 1         | -            |
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # Stranded PyRanges object has 10 rows and 6 columns from 1 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("aorta2.bed")

    return pr.read_bed(full_path)


def bw():
    full_path = get_example_path("bw.bw")

    return pr.read_bigwig(full_path)


def chipseq():
    """
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # | Chromosome   | Start     | End       | Name       | Score     | Strand       |
    >>> # | (category)   | (int64)   | (int64)   | (object)   | (int64)   | (category)   |
    >>> # |--------------+-----------+-----------+------------+-----------+--------------|
    >>> # | chr1         | 212609534 | 212609559 | U0         | 0         | +            |
    >>> # | chr1         | 169887529 | 169887554 | U0         | 0         | +            |
    >>> # | chr1         | 216711011 | 216711036 | U0         | 0         | +            |
    >>> # | chr1         | 144227079 | 144227104 | U0         | 0         | +            |
    >>> # | ...          | ...       | ...       | ...        | ...       | ...          |
    >>> # | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
    >>> # | chrY         | 13517892  | 13517917  | U0         | 0         | -            |
    >>> # | chrY         | 8010951   | 8010976   | U0         | 0         | -            |
    >>> # | chrY         | 7405376   | 7405401   | U0         | 0         | -            |
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # Stranded PyRanges object has 10,000 rows and 6 columns from 24 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("chipseq.bed")

    return pr.read_bed(full_path)


def chipseq_background():
    """
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # | Chromosome   | Start     | End       | Name       | Score     | Strand       |
    >>> # | (category)   | (int64)   | (int64)   | (object)   | (int64)   | (category)   |
    >>> # |--------------+-----------+-----------+------------+-----------+--------------|
    >>> # | chr1         | 39036822  | 39036847  | U0         | 0         | +            |
    >>> # | chr1         | 224145989 | 224146014 | U0         | 0         | +            |
    >>> # | chr1         | 167802964 | 167802989 | U0         | 0         | +            |
    >>> # | chr1         | 69101066  | 69101091  | U0         | 0         | +            |
    >>> # | ...          | ...       | ...       | ...        | ...       | ...          |
    >>> # | chrY         | 11936866  | 11936891  | U0         | 0         | -            |
    >>> # | chrY         | 10629111  | 10629136  | U0         | 0         | -            |
    >>> # | chrY         | 10632456  | 10632481  | U0         | 0         | -            |
    >>> # | chrY         | 11918814  | 11918839  | U0         | 0         | -            |
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # Stranded PyRanges object has 10,000 rows and 6 columns from 25 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("chipseq_background.bed")

    return pr.read_bed(full_path)


def chromsizes():
    """
    >>> # +--------------+-----------+-----------+
    >>> # | Chromosome   | Start     | End       |
    >>> # | (category)   | (int64)   | (int64)   |
    >>> # |--------------+-----------+-----------|
    >>> # | chr1         | 0         | 249250621 |
    >>> # | chr2         | 0         | 243199373 |
    >>> # | chr3         | 0         | 198022430 |
    >>> # | chr4         | 0         | 191154276 |
    >>> # | ...          | ...       | ...       |
    >>> # | chrY         | 0         | 59373566  |
    >>> # | chrX         | 0         | 155270560 |
    >>> # | chrM         | 0         | 16571     |
    >>> # | chr22        | 0         | 51304566  |
    >>> # +--------------+-----------+-----------+
    >>> # Unstranded PyRanges object has 25 rows and 3 columns from 25 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome.
    """

    full_path = get_example_path("chromsizes.bed")

    return pr.read_bed(full_path)


def control_bam():
    """
    >>> # +--------------+-----------+-----------+--------------+------------+
    >>> # | Chromosome   | Start     | End       | Strand       | Flag       |
    >>> # | (category)   | (int64)   | (int64)   | (category)   | (uint16)   |
    >>> # |--------------+-----------+-----------+--------------+------------|
    >>> # | chr1         | 887771    | 887796    | +            | 16         |
    >>> # | chr1         | 994660    | 994685    | +            | 16         |
    >>> # | chr1         | 1770383   | 1770408   | +            | 16         |
    >>> # | chr1         | 1995141   | 1995166   | +            | 16         |
    >>> # | ...          | ...       | ...       | ...          | ...        |
    >>> # | chrY         | 57402214  | 57402239  | +            | 16         |
    >>> # | chrY         | 10643526  | 10643551  | -            | 0          |
    >>> # | chrY         | 11776321  | 11776346  | -            | 0          |
    >>> # | chrY         | 20557165  | 20557190  | -            | 0          |
    >>> # +--------------+-----------+-----------+--------------+------------+
    >>> # Stranded PyRanges object has 10,000 rows and 5 columns from 25 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("control.bam")

    return pr.read_bam(full_path)


def cpg():
    """
    >>> # +--------------+-----------+-----------+-----------+
    >>> # | Chromosome   | Start     | End       | CpG       |
    >>> # | (category)   | (int64)   | (int64)   | (int64)   |
    >>> # |--------------+-----------+-----------+-----------|
    >>> # | chrX         | 64181     | 64793     | 62        |
    >>> # | chrX         | 69133     | 70029     | 100       |
    >>> # | chrX         | 148685    | 149461    | 85        |
    >>> # | chrX         | 166504    | 167721    | 96        |
    >>> # | ...          | ...       | ...       | ...       |
    >>> # | chrY         | 28555535  | 28555932  | 32        |
    >>> # | chrY         | 28773315  | 28773544  | 25        |
    >>> # | chrY         | 59213794  | 59214183  | 36        |
    >>> # | chrY         | 59349266  | 59349574  | 29        |
    >>> # +--------------+-----------+-----------+-----------+
    >>> # Unstranded PyRanges object has 1,077 rows and 4 columns from 2 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome.
    """

    full_path = get_example_path("cpg.bed")

    df = pd.read_csv(full_path, sep="\t", header=None, names="Chromosome Start End CpG".split())

    return pr.PyRanges(df)


def ensembl_gtf():
    """
    >>> # +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------+
    >>> # | Chromosome   | Source     | Feature      | Start     | End       | Score      | Strand       | Frame      | gene_biotype                       | +19   |
    >>> # | (category)   | (object)   | (category)   | (int64)   | (int64)   | (object)   | (category)   | (object)   | (object)                           | ...   |
    >>> # |--------------+------------+--------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------|
    >>> # | 1            | havana     | gene         | 11868     | 14409     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
    >>> # | 1            | havana     | transcript   | 11868     | 14409     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
    >>> # | 1            | havana     | exon         | 11868     | 12227     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
    >>> # | 1            | havana     | exon         | 12612     | 12721     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
    >>> # | ...          | ...        | ...          | ...       | ...       | ...        | ...          | ...        | ...                                | ...   |
    >>> # | 1            | havana     | gene         | 1173055   | 1179555   | .          | -            | .          | lncRNA                             | ...   |
    >>> # | 1            | havana     | transcript   | 1173055   | 1179555   | .          | -            | .          | lncRNA                             | ...   |
    >>> # | 1            | havana     | exon         | 1179364   | 1179555   | .          | -            | .          | lncRNA                             | ...   |
    >>> # | 1            | havana     | exon         | 1173055   | 1176396   | .          | -            | .          | lncRNA                             | ...   |
    >>> # +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------+
    >>> # Stranded PyRanges object has 2,446 rows and 28 columns from 1 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    >>> # 19 hidden columns: gene_id, gene_name, gene_source, gene_version, tag, transcript_biotype, transcript_id, transcript_name, transcript_source, transcript_support_level, ... (+ 9 more.)
    """

    full_path = get_example_path("ensembl_human.gtf.gz")

    return pr.read_gtf(full_path)


def exons():
    """
    >>> # +--------------+-----------+-----------+----------------------------------------+-----------+--------------+
    >>> # | Chromosome   | Start     | End       | Name                                   | Score     | Strand       |
    >>> # | (category)   | (int64)   | (int64)   | (object)                               | (int64)   | (category)   |
    >>> # |--------------+-----------+-----------+----------------------------------------+-----------+--------------|
    >>> # | chrX         | 135721701 | 135721963 | NR_038462_exon_0_0_chrX_135721702_f    | 0         | +            |
    >>> # | chrX         | 135574120 | 135574598 | NM_001727_exon_2_0_chrX_135574121_f    | 0         | +            |
    >>> # | chrX         | 47868945  | 47869126  | NM_205856_exon_4_0_chrX_47868946_f     | 0         | +            |
    >>> # | chrX         | 77294333  | 77294480  | NM_000052_exon_17_0_chrX_77294334_f    | 0         | +            |
    >>> # | ...          | ...       | ...       | ...                                    | ...       | ...          |
    >>> # | chrY         | 15409586  | 15409728  | NR_047633_exon_3_0_chrY_15409587_r     | 0         | -            |
    >>> # | chrY         | 15478146  | 15478273  | NR_047634_exon_18_0_chrY_15478147_r    | 0         | -            |
    >>> # | chrY         | 15360258  | 15361762  | NR_047601_exon_0_0_chrY_15360259_r     | 0         | -            |
    >>> # | chrY         | 15467254  | 15467278  | NM_001258270_exon_13_0_chrY_15467255_r | 0         | -            |
    >>> # +--------------+-----------+-----------+----------------------------------------+-----------+--------------+
    >>> # Stranded PyRanges object has 1,000 rows and 6 columns from 2 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("exons.bed")

    return pr.read_bed(full_path)


def f1():
    """
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # | Chromosome   |     Start |       End | Name       |     Score | Strand       |
    >>> # | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
    >>> # |--------------+-----------+-----------+------------+-----------+--------------|
    >>> # | chr1         |         3 |         6 | interval1  |         0 | +            |
    >>> # | chr1         |         8 |         9 | interval3  |         0 | +            |
    >>> # | chr1         |         5 |         7 | interval2  |         0 | -            |
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("f1.bed")

    return pr.read_bed(full_path)


def f2():
    """
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # | Chromosome   |     Start |       End | Name       |     Score | Strand       |
    >>> # | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
    >>> # |--------------+-----------+-----------+------------+-----------+--------------|
    >>> # | chr1         |         1 |         2 | a          |         0 | +            |
    >>> # | chr1         |         6 |         7 | b          |         0 | -            |
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("f2.bed")

    return pr.read_bed(full_path)


def gencode_gtf():
    """
    >>> # +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+-------------------+-------+
    >>> # | Chromosome   | Source     | Feature      | Start     | End       | Score      | Strand       | Frame      | gene_id           | +15   |
    >>> # | (category)   | (object)   | (category)   | (int64)   | (int64)   | (object)   | (category)   | (object)   | (object)          | ...   |
    >>> # |--------------+------------+--------------+-----------+-----------+------------+--------------+------------+-------------------+-------|
    >>> # | chr1         | HAVANA     | gene         | 11868     | 14409     | .          | +            | .          | ENSG00000223972.5 | ...   |
    >>> # | chr1         | HAVANA     | transcript   | 11868     | 14409     | .          | +            | .          | ENSG00000223972.5 | ...   |
    >>> # | chr1         | HAVANA     | exon         | 11868     | 12227     | .          | +            | .          | ENSG00000223972.5 | ...   |
    >>> # | chr1         | HAVANA     | exon         | 12612     | 12721     | .          | +            | .          | ENSG00000223972.5 | ...   |
    >>> # | ...          | ...        | ...          | ...       | ...       | ...        | ...          | ...        | ...               | ...   |
    >>> # | chr1         | HAVANA     | exon         | 1430549   | 1430662   | .          | -            | .          | ENSG00000225285.1 | ...   |
    >>> # | chr1         | HAVANA     | transcript   | 1430663   | 1434520   | .          | -            | .          | ENSG00000225285.1 | ...   |
    >>> # | chr1         | HAVANA     | exon         | 1434177   | 1434520   | .          | -            | .          | ENSG00000225285.1 | ...   |
    >>> # | chr1         | HAVANA     | exon         | 1430663   | 1430954   | .          | -            | .          | ENSG00000225285.1 | ...   |
    >>> # +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+-------------------+-------+
    >>> # Stranded PyRanges object has 4,995 rows and 24 columns from 1 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    >>> # 15 hidden columns: gene_type, gene_name, level, havana_gene, transcript_id, transcript_type, transcript_name, transcript_support_level, tag, ... (+ 6 more.)
    """

    full_path = get_example_path("gencode_human.gtf.gz")

    return pr.read_gtf(full_path)


def ucsc_bed():
    """
    >>> # +--------------+-----------+-----------+------------+------------+-----------------+--------------+---------------+-------------------+
    >>> # | Chromosome   | Start     | End       | Feature    | gene_id    | transcript_id   | Strand       | exon_number   | transcript_name   |
    >>> # | (category)   | (int64)   | (int64)   | (object)   | (object)   | (float64)       | (category)   | (float64)     | (object)          |
    >>> # |--------------+-----------+-----------+------------+------------+-----------------+--------------+---------------+-------------------|
    >>> # | chr1         | 12776117  | 12788726  | gene       | AADACL3    | nan             | +            | nan           | nan               |
    >>> # | chr1         | 169075927 | 169101957 | gene       | ATP1B1     | nan             | +            | nan           | nan               |
    >>> # | chr1         | 6845383   | 7829766   | gene       | CAMTA1     | nan             | +            | nan           | nan               |
    >>> # | chr1         | 20915589  | 20945396  | gene       | CDA        | nan             | +            | nan           | nan               |
    >>> # | ...          | ...       | ...       | ...        | ...        | ...             | ...          | ...           | ...               |
    >>> # | chrX         | 152661096 | 152663330 | exon       | PNMA6E     | 260.0           | -            | 0.0           | NM_001351293      |
    >>> # | chrX         | 152661096 | 152666808 | transcript | PNMA6E     | 260.0           | -            | nan           | NM_001351293      |
    >>> # | chrX         | 152664164 | 152664378 | exon       | PNMA6E     | 260.0           | -            | 1.0           | NM_001351293      |
    >>> # | chrX         | 152666701 | 152666808 | exon       | PNMA6E     | 260.0           | -            | 2.0           | NM_001351293      |
    >>> # +--------------+-----------+-----------+------------+------------+-----------------+--------------+---------------+-------------------+
    >>> # Stranded PyRanges object has 5,519 rows and 9 columns from 30 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("ucsc_human.bed.gz")

    names = "Chromosome Start End Feature gene_id transcript_id Strand exon_number transcript_name".split()
    df = pd.read_csv(full_path, sep="\t", names=names)

    return pr.PyRanges(df)
