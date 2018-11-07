import pytest
from tests.helpers import assert_dfs_equal, string_to_pyrange

from pyranges.pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO

@pytest.fixture
def expected_result_set_intersection():

    c = """Chromosome Start End
chr1	226987603	226987617
chr8	38747236	38747251
chr15	26105515	26105518"""

    return string_to_pyrange(c)


def test_advanced_intersection(cs, bg, expected_result_set_intersection):

    result = cs.set_intersection(bg, strandedness=False)

    print(result)
    print(expected_result_set_intersection)

    assert_dfs_equal(result, expected_result_set_intersection)


@pytest.fixture
def expected_result_set_intersection_same_strand():

    c = """Chromosome Start End Strand
chr15	26105515	26105518 +"""

    return string_to_pyrange(c)

def test_advanced_intersection_same_strand(cs, bg, expected_result_set_intersection_same_strand):

    result = cs.set_intersection(bg, strandedness="same")

    assert_dfs_equal(result, expected_result_set_intersection_same_strand)


# @pytest.fixture
# def expected_result_intersection_same_strand():

#     c = """chr15	26105515	26105518	U0	0	+"""

#     return pd.read_table(StringIO(c), header=None, names="Chromosome Start End Name Score Strand".split())

# def test_advanced_intersection_stranded(cs, bg, expected_result_intersection_same_strand):

#     result = cs.intersection(bg, strandedness="same")

#     print(result.df)
#     print(expected_result_intersection_same_strand)

#     assert_df_equal(expected_result_intersection_same_strand, result.df)


# @pytest.fixture
# def expected_result_intersection_opposite_strand():

#     c = """chr1	226987603	226987617	U0	0	+
# chr8	38747236	38747251	U0	0	-"""

#     return pd.read_table(StringIO(c), header=None, names="Chromosome Start End Name Score Strand".split())


# def test_advanced_intersection_opposite_stranded(cs, bg, expected_result_intersection_opposite_strand):

#     result = cs.intersection(bg, strandedness="opposite")

#     print("result")
#     print(result.df)
#     print("expected result")
#     print(expected_result_intersection_opposite_strand)

#     assert_df_equal(expected_result_intersection_opposite_strand, result.df)
