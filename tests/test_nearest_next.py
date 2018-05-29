import pytest
from tests.helpers import assert_df_equal

from pyranges.pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO


@pytest.fixture
def expected_result_unstranded():

 c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
0 chr1 3 6 h 0 + 6 7 f 0 - 0
1 chr1 5 7 h 0 - 6 7 f 0 - 0"""

 return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_unstranded(f1, f2, expected_result_unstranded):

    result = f1.nearest(f2, strandedness=False, suffix="_b", how="next")

    print(result.df.to_csv(sep=" "))

    assert assert_df_equal(result.df, expected_result_unstranded.df)


@pytest.fixture
def expected_result_opposite_stranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
0 chr1 3 6 h 0 + 6 7 f 0 - 0"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_opposite_stranded(f1, f2, expected_result_opposite_stranded):

    result = f1.nearest(f2, strandedness="opposite", suffix="_b", how="next")

    print(result.df.to_csv(sep=" "))

    assert assert_df_equal(result.df, expected_result_opposite_stranded.df)


@pytest.fixture
def expected_result_same_stranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
0 chr1 5 7 h 0 - 6 7 f 0 - 0"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_same_stranded(f1, f2, expected_result_same_stranded):

    result = f1.nearest(f2, strandedness="same", suffix="_b", how="next")

    print(result.df.to_csv(sep=" "))

    assert assert_df_equal(result.df, expected_result_same_stranded.df)
