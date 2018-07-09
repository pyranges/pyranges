

import pytest
from tests.helpers import assert_df_equal

from pyranges.pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO


def test_tss(f1):

    result = f1.tss(slack=5)

    print(f1)
    print(result)

    assert list(result.Start) == [0, 2, 3]
    assert list(result.End) == [9, 13, 14]


def test_tes(f2):

    print(f2)

    result = f2.tes(slack=500)

    print(result)

    assert list(result.Start) == [0, 0]
    assert list(result.End) == [503, 508]
