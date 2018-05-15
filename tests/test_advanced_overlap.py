import pytest

from tests.helpers import assert_df_equal
from pyranges.pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO

@pytest.fixture
def cs():

    cs = pr.load_dataset("chipseq")
    return cs



@pytest.fixture
def bg():

    bg = pr.load_dataset("chipseq_background")
    return bg


@pytest.fixture
def expected_result_regular_overlap():

    c = """chr1	226987592	226987617	U0	0	+
chr8	38747226	38747251	U0	0	-
chr15	26105515	26105540	U0	0	+"""

    return pd.read_table(StringIO(c), header=None, names="Chromosome Start End Name Score Strand".split())


def test_advanced_overlap(cs, bg, expected_result_regular_overlap):

    result = cs.overlap(bg)

    print(result.df)
    print(expected_result_regular_overlap)

    assert assert_df_equal(result.df, expected_result_regular_overlap)

@pytest.fixture
def expected_result_overlap_same_strand():

    c = """chr15	26105515	26105540	U0	0	+"""

    return pd.read_table(StringIO(c), header=None, names="Chromosome Start End Name Score Strand".split())

def test_advanced_overlap_stranded(cs, bg, expected_result_overlap_same_strand):

    result = cs.overlap(bg, strandedness="same")

    print(result.df)
    print(expected_result_overlap_same_strand)

    assert assert_df_equal(result.df, expected_result_overlap_same_strand)


@pytest.fixture
def expected_result_overlap_opposite_strand():

    c = """chr1	226987592	226987617	U0	0	+
chr8	38747226	38747251	U0	0	-"""

    return pd.read_table(StringIO(c), header=None, names="Chromosome Start End Name Score Strand".split())


def test_advanced_overlap_opposite_stranded(cs, bg, expected_result_overlap_opposite_strand):

    result = cs.overlap(bg, strandedness="opposite")

    print(result.df)
    print(expected_result_overlap_opposite_strand)

    assert assert_df_equal(result.df, expected_result_overlap_opposite_strand)
