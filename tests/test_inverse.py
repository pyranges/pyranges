
import pytest

import pandas as pd
import numpy as np

from pyrle import Rle



@pytest.fixture
def rle():

    r = pd.Series([0], dtype=np.int16)
    r1 = Rle(r, r)

    return r1

# @pytest.fixture
# def rle2():

#     r = pd.Series([2, 1], dtype=np.int16)
#     # v = pd.Series([0, 1], dtype=np.int16)
#     r2 = Rle(r, v)

#     return r2

def test_hypothesis_counterexample(rle):

    result = rle + rle
    result = result - rle

    expected_runs = []
    expected_values = []

    print(result.runs)
    print(result.values)

    assert list(result.runs) == expected_runs
    assert list(result.values) == expected_values

    # result2 = rle2 + rle

    # print(result2.runs)
    # print(result2.values)

    # assert list(result2.runs) == expected_runs
    # assert list(result2.values) == expected_values
