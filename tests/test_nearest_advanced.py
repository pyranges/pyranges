
import pytest
from tests.helpers import assert_df_equal

from pyranges.pyranges import PyRanges


import pandas as pd

from io import StringIO

@pytest.fixture
def expected_result_unstranded():
    pass


# def test_advanced_nearest(cs, bg):

#     result = cs.nearest(bg, suffix="_b")

#     print(result)

#     assert 0
