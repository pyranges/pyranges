

import pytest
from tests.helpers import assert_df_equal

from pyranges.pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO


def test_tssify(f1):

    result = f1.tssify(slack=5)

    print(f1)
    print(result)

    assert list(result.Start) == [0, 2, 3]
    assert list(result.End) == [9, 13, 14]


def test_tesify(f2):

    print(f2)

    result = f2.tesify(slack=500)

    print(result)

    assert list(result.Start) == [0, 0]
    assert list(result.End) == [503, 508]
