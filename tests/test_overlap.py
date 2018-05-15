
from tests.helpers import assert_df_equal
import pytest

from pyranges.pyranges import PyRanges

import pandas as pd

from io import StringIO


@pytest.fixture
def simple_gr1():

    c = """Chromosome Start End Strand Score
chr1 3 6 + 5
chr1 5 7 - 7
chr1 8 9 + 1"""

    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return PyRanges(df)



@pytest.fixture
def simple_gr2():

    c = """Chromosome Start End Strand Score
chr1 1 2 + 1
chr1 6 7 - 2"""
    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return PyRanges(df)


@pytest.fixture
def expected_result_subtract_simple_granges():

    c = """Chromosome Start End Strand Score
chr1	3	6	+ 5
chr1	8	9	+ 1"""
    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return PyRanges(df)



def test_overlap_invert_simple_granges(simple_gr1, simple_gr2, expected_result_subtract_simple_granges):

    result = simple_gr1.overlap(simple_gr2, invert=True, strandedness=False)

    print(result)
    print(expected_result_subtract_simple_granges)


    print(result.df.dtypes)
    print(expected_result_subtract_simple_granges.df.dtypes)

    assert assert_df_equal(result.df, expected_result_subtract_simple_granges.df)



@pytest.fixture
def expected_result_overlap_same_strand_simple_granges():

    c = """Chromosome Start End Strand Score
chr1 5 7 - 7"""
    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return PyRanges(df)


def test_overlap_same_strand_simple_granges(simple_gr1, simple_gr2, expected_result_overlap_same_strand_simple_granges):

    result = simple_gr1.overlap(simple_gr2, strandedness="same")

    print(result)

    assert assert_df_equal(result.df, expected_result_overlap_same_strand_simple_granges.df)



def test_overlap_opposite_strand_simple_granges(simple_gr1, simple_gr2):

    result = simple_gr1.overlap(simple_gr2, strandedness="opposite")

    print(result)

    assert result.df.empty


def test_default_overlap_simple_granges(simple_gr1, simple_gr2, expected_result_overlap_same_strand_simple_granges):

    print(simple_gr1)
    print(simple_gr2)

    result = simple_gr1.overlap(simple_gr2)

    print(result.df.dtypes)
    print(expected_result_overlap_same_strand_simple_granges.df.dtypes)

    print(result.df.index)
    print(expected_result_overlap_same_strand_simple_granges.df.index)

    assert assert_df_equal(result.df, expected_result_overlap_same_strand_simple_granges.df)
