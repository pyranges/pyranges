
import pytest

from pyranges.pyranges import PyRanges

from tests.helpers import assert_df_equal
import pandas as pd
import numpy as np

from io import StringIO

def test_coverage_score(cs):

    result = cs.coverage(value_col="Score")

    assert result["chr1"].runs == np.array([247134924], dtype=np.long)


def test_coverage(cs):

    result = cs.coverage()

    print(list(result["chrY"].runs))

    assert list(result["chrY"].runs) == [7046809, 25, 147506, 25, 211011, 25, 58043, 25, 238514, 25, 59018, 25, 249900, 25, 305797, 25, 3625972, 25, 987578, 25, 286216, 25, 301253, 25, 1256136, 25, 450157, 25, 323762, 25, 497195, 25, 450230, 25, 5063659, 25, 148456, 25, 43524, 25, 159470, 25, 143271, 25, 156610, 25]
