import pytest

import pandas as pd

from io import StringIO

from pyranges import PyRanges

from tests.helpers import assert_df_equal

@pytest.fixture
def expected_result_previous_bed_opposite_stranded(names):

    c = """chr1 3 6 h 0 + 6 7 f 0 - 0
chr1 8 9 h 0 + 6 7 f 0 - 1
chr1 5 7 h 0 - 1 2 f 0 + 3"""

    df = pd.read_table(StringIO(c), sep=" ", header=None, names="Chromosome  Start  End  Name Score Strand Start_b  End_b  Name_b Score_b Strand_b Distance".split())
    print(df)

    return PyRanges(df)

def test_nearest_previous_bed_opposite_stranded(f1, f2, expected_result_previous_bed_opposite_stranded):

    print(f1)
    print(f2)
    print(expected_result_previous_bed_opposite_stranded)
    result = f1.nearest(f2, strandedness="opposite", suffix="_b", how="previous")

    print(result.df.to_csv(sep=" "))

    assert assert_df_equal(result.df, expected_result_previous_bed_opposite_stranded.df)


@pytest.fixture
def expected_result_previous_bed_same_stranded(names):

    c = """0 chr1 3 6 h 0 + 1 2 f 0 + 1
1 chr1 5 7 h 0 - 6 7 f 0 - 0
2 chr1 8 9 h 0 + 1 2 f 0 + 6"""

    df = pd.read_table(StringIO(c), sep=" ", header=None, names="Chromosome  Start  End  Name Score Strand Start_b  End_b  Name_b Score_b Strand_b Distance".split())
    print(df)

    return PyRanges(df)


def test_nearest_previous_bed_same_stranded(f1, f2, expected_result_previous_bed_same_stranded):

    result = f1.nearest(f2, strandedness="same", suffix="_b", how="previous")

    assert assert_df_equal(result.df, expected_result_previous_bed_same_stranded.df)


@pytest.fixture
def expected_result_previous_bed_unstranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
0 chr1 3 6 h 0 + 6 7 f 0 - 0
1 chr1 5 7 h 0 - 6 7 f 0 - 0
2 chr1 8 9 h 0 + 6 7 f 0 - 1"""

    df = pd.read_table(StringIO(c), sep=" ", header=0)
    print(df)

    return PyRanges(df)

def test_nearest_previous_bed_unstranded(f1, f2, expected_result_previous_bed_unstranded):

    result = f1.nearest(f2, strandedness=None, suffix="_b", how="previous")

    assert assert_df_equal(result.df, expected_result_previous_bed_unstranded.df)
