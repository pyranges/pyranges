import pytest
from tests.helpers import assert_df_equal

from pyranges.pyranges import PyRanges
from pyranges.methods import both_indexes
import pyranges as pr

import pandas as pd
from collections import OrderedDict

from io import StringIO

@pytest.fixture
def expected_result_both_indexes():

    return OrderedDict([("chr1", [0, 0])]), OrderedDict([("chr1", [0, 1])])


@pytest.fixture
def c1(names):

    c = """chr1	1	2	HWI-ST216_313:3:1203:10227:6568	1	+
chr1	2	3	HWI-ST216_313:3:1203:10227:6568	1	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", names=names, header=None))

@pytest.fixture
def c2(names):

    c = """chr1	1	3	HWI-ST216_313:3:1203:10227:6568	1	+
chr1	1	2	HWI-ST216_313:3:1203:10227:6568	1	-"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", names=names, header=None))




def test_both_indexes(c1, c2, expected_result_both_indexes):

    sidx, oidx = both_indexes(c1, c2, strandedness=False, how="first")

    print(sidx)
    print(oidx)

    assert result == expected_result_both_indexes


# do afterwards
# def test_intersect_bed_opposite_strand_containment(c1, c2, expected_result_both_indexes):
#     result = both_indexes(c1, c2, strandedness=False, how="first")
#     assert result == expected_result_both_indexes
