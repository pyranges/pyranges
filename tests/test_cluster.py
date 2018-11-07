import pytest
from tests.helpers import assert_dfs_equal

from pyranges.pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO


@pytest.fixture
def expected_result_cluster_simple_granges():

    c = """Chromosome Start End Strand
chr1 3 6 +
chr1 5 7 -
chr1 8 9 +"""

    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return PyRanges(df)

@pytest.fixture
def expected_result_cluster_simple_granges_no_strand():

    c = """Chromosome Start End
chr1 3 7
chr1 8 9"""

    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return PyRanges(df)



# @pytest.fixture

#     c = """Chromosome Start End Strand Score
# chr1    6   7   - 7"""

#     df = pd.read_table(StringIO(c), sep="\s+", header=0)
#     return PyRanges(df)


def test_cluster_simple_granges(f1, expected_result_cluster_simple_granges):

    result = f1.cluster(strand=True)

    print("result\n", result)

    assert_dfs_equal(result, expected_result_cluster_simple_granges)


def test_cluster_simple_granges_no_strand(f1, expected_result_cluster_simple_granges_no_strand):

    result = f1.cluster(strand=None)

    print("result\n", result)

    assert_dfs_equal(result, expected_result_cluster_simple_granges_no_strand)
