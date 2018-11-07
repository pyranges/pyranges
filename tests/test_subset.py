import pytest

from pyranges.pyranges import PyRanges

from tests.helpers import assert_dfs_equal
import pandas as pd

from io import StringIO


# @pytest.fixture
# def expected_result_subset():

#     c = """Chromosome Start End Strand
# chr1 3 6 +
# chr1 5 7 -
# chr1 8 9 +"""

#     df = pd.read_table(StringIO(c), sep="\s+", header=0)
#     return PyRanges(df)

# @pytest.fixture
# def expected_result_cluster_simple_granges_no_strand():

#     c = """Chromosome Start End
# chr1 3 7
# chr1 8 9"""

#     df = pd.read_table(StringIO(c), sep="\s+", header=0)
#     return PyRanges(df)



# @pytest.fixture

#     c = """Chromosome Start End Strand Score
# chr1    6   7   - 7"""

#     df = pd.read_table(StringIO(c), sep="\s+", header=0)
#     return PyRanges(df)


def test_subset_chr(cs):

    # result = pd.concat(list(cs["chrY"].values()))
    cs["chrY", 100:100000]

    print(result)

    assert 0
    # assert_dfs_equal(result)
