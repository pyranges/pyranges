import pyranges as pr

# def test_read_bam():

#     pr.read_bam("tests/test_data/test_sorted.bam")


def test_read_gtf():

    gr = pr.read_gtf("tests/test_data/ensembl.gtf", full=True)
    assert len(gr.columns) == 28

    # assert list(gr.df.columns[:4]) == "Chromosome Start End Strand".split()


def test_read_gff3():

    gr = pr.read_gff3("tests/test_data/gencode.gff3")

    assert len(gr.columns) == 26
    # assert list(gr.df.columns[:4]) == "Chromosome Start End Strand".split()


def test_read_bed():

    pr.read_bed("pyranges/example_data/chipseq.bed")
