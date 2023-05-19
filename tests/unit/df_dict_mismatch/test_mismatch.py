import pandas as pd

import pyranges as pr


def test_unstranded_but_has_chrom_key():
    df = pd.DataFrame({"Chromosome": "chr1", "Start": 5, "End": 10}, index=[0])
    dfs = {("chr1", "+"): df}

    gr = pr.PyRanges(dfs)

    assert not gr.stranded


def test_has_bad_strand_and_strand_key():
    df = pd.DataFrame({"Chromosome": "chr1", "Start": 5, "End": 10, "Strand": "."}, index=[0])

    dfs = {("chr1", "+"): df}

    gr = pr.PyRanges(dfs)

    assert not gr.stranded


def test_has_strand_but_is_not_stranded():
    df = pd.DataFrame({"Chromosome": "chr1", "Start": 5, "End": 10, "Strand": "+"}, index=[0])

    dfs = {("chr1"): df}

    gr = pr.PyRanges(dfs)

    print(gr.dfs)

    assert gr.stranded
