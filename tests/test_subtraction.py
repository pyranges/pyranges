import pytest

from tests.helpers import assert_df_equal
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
def expected_result_inverse_subtraction_simple_granges():

    c = """Chromosome Start End Strand Score
chr1	3	6	+ 5
chr1    5   6   - 7
chr1	8	9	+ 1"""

    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return PyRanges(df)



def test_intersect_invert_simple_granges(simple_gr1, simple_gr2, expected_result_inverse_subtraction_simple_granges):

    result = simple_gr1.subtraction(simple_gr2, strandedness=False)

    assert assert_df_equal(result.df, expected_result_inverse_subtraction_simple_granges.df)


@pytest.fixture
def expected_result_same_strand_inverse_subtraction_simple_granges():

    c = """Chromosome Start End Strand Score
chr1	3	6	+ 5
chr1    5   6   - 7
chr1	8	9	+ 1"""

    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return PyRanges(df)


def test_intersect_invert_same_strand_simple_granges(simple_gr1, simple_gr2, expected_result_same_strand_inverse_subtraction_simple_granges):

    result = simple_gr1.subtraction(simple_gr2, strandedness="same")

    print(result)

    assert assert_df_equal(result.df, expected_result_same_strand_inverse_subtraction_simple_granges.df)



@pytest.fixture
def expected_result_opposite_strand_inverse_subtraction_simple_granges():

    c = """Chromosome Start End Strand Score
chr1	3	6	+ 5
chr1    5   7   - 7
chr1	8	9	+ 1"""

    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return PyRanges(df)


def test_intersect_invert_opposite_strand_simple_granges(simple_gr1, simple_gr2, expected_result_opposite_strand_inverse_subtraction_simple_granges):

    result = simple_gr1.subtraction(simple_gr2, strandedness="opposite")

    print(result)

    assert assert_df_equal(result.df, expected_result_opposite_strand_inverse_subtraction_simple_granges.df)
