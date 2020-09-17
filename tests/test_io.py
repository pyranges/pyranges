import pyranges as pr


def test_read_gtf():

    gr = pr.read_gtf("tests/test_data/ensembl.gtf", full=True)

    assert len(gr.columns) == 26

    df = gr.df
    transcript = df.iloc[1]
    assert transcript['tag'] == 'basic'

    exon = df[df['exon_id'] == 'ENSE00003812156'].iloc[0]
    assert exon['tag'] == 'basic'

    gr = pr.read_gtf("tests/test_data/ensembl.gtf",
                     full=True, duplicate_attr=True)
    print(gr.columns)
    assert len(gr.columns) == 26

    df = gr.df
    transcript = df.iloc[1]
    assert transcript['tag'] == 'basic'

    exon = df[df['exon_id'] == 'ENSE00003812156'].iloc[0]
    assert exon['tag'] == 'CCDS,basic'
    # assert list(gr.df.columns[:4]) == "Chromosome Start End Strand".split()


def test_read_gff3():

    gr = pr.read_gff3("tests/test_data/gencode.gff3")

    assert len(gr.columns) == 26
    # assert list(gr.df.columns[:4]) == "Chromosome Start End Strand".split()


def test_read_bed():

    pr.read_bed("pyranges/example_data/chipseq.bed")
