from pyranges.pyranges import PyRanges

from tests.helpers import assert_df_equal
import pytest

from io import StringIO

import pandas as pd


@pytest.fixture
def background():
    c = """Chromosome Start End Name Score Strand
chr1	226987603	226987628	U0	0	-
chr8	38747236	38747261	U0	0	+
chr15	26105493	26105518	U0	0	+"""


    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    print(df)
    return PyRanges(df)

@pytest.fixture
def chip():

    c = """Chromosome Start End Name Score Strand
chr1	226987592	226987617	U0	0	+
chr8	38747226	38747251	U0	0	-
chr15	26105515	26105540	U0	0	+"""


    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    # print(df)
    return PyRanges(df)



@pytest.fixture
def expected_result_regular_intersection():

    c = """chr1	226987603	226987617	U0	0	+
chr8	38747236	38747251	U0	0	-
chr15	26105515	26105518	U0	0	+"""

    return PyRanges(pd.read_table(StringIO(c), header=None, names="Chromosome Start End Name Score Strand".split()))


def test_intersect_granges(chip, background, expected_result_regular_intersection):

    print("chip " * 3)
    print(chip)

    print("background " * 3)
    print(background)

    print("expected_result " * 3)
    print(expected_result_regular_intersection)

    result = chip.intersection(background, strandedness=False)

    print(result)

    assert assert_df_equal(result.df, expected_result_regular_intersection.df)
