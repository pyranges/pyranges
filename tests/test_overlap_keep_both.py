import pytest

from tests.helpers import assert_df_equal
from pyranges.pyranges import PyRanges
from pyranges.methods import _overlap_write_both

import pandas as pd

from io import StringIO


@pytest.fixture
def simple_gr1():

    c = """Chromosome Start End Score Strand
chr1 3 6 5 +
chr1 5 7 7 -
chr1 8 9 1 +"""

    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return PyRanges(df)



@pytest.fixture
def simple_gr2():

    c = """Chromosome Start End Score Strand
chr1 1 2 1 +
chr1 6 7 2 -"""
    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return PyRanges(df)


@pytest.fixture
def expected_result_subtract_simple_granges():

    c = """Chromosome Start End Strand Score
chr1	3	6	+ 5
chr1	8	9	+ 1"""
    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return PyRanges(df)




@pytest.fixture
def expected_result_overlap_same_strand_simple_granges():

    c = """Chromosome Start End Score Strand Start_b End_b Score_b Strand_b
chr1	5	7	7	-	6	7	2	-"""
    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return PyRanges(df)


def test_overlap_same_strand_simple_granges(simple_gr1, simple_gr2, expected_result_overlap_same_strand_simple_granges):

    print("gr1")
    print(simple_gr1)

    print("gr2")
    print(simple_gr2)

    print("expected")
    print(expected_result_overlap_same_strand_simple_granges)
    result = simple_gr1.join(simple_gr2, strandedness="same")

    print("actual")
    print(result)

    assert assert_df_equal(result.df, expected_result_overlap_same_strand_simple_granges.df)



def test_overlap_opposite_strand_simple_granges(simple_gr1, simple_gr2):

    result = simple_gr1.join(simple_gr2, strandedness="opposite")

    print(result)

    assert result.df.empty


def test_default_overlap_simple_granges(simple_gr1, simple_gr2, expected_result_overlap_same_strand_simple_granges):

    print(simple_gr1)
    print(simple_gr2)

    result = simple_gr1.join(simple_gr2)

    assert assert_df_equal(result.df, expected_result_overlap_same_strand_simple_granges.df)
