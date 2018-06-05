

import pytest

import pandas as pd
import numpy as np

from pyrle import Rle

@pytest.fixture
def weird_rle():

    return Rle([10, 20, 30, 40], [1, 2, 3, 4])



@pytest.fixture
def weird_rle2():

    return Rle([1], [1])


@pytest.fixture
def expected_result_weird():

    return Rle([1, 99], [1, 0])


def test_weird(weird_rle, weird_rle2, expected_result_weird):

    result = weird_rle2 * weird_rle
    result2 = weird_rle * weird_rle2

    print(result)

    assert result == expected_result_weird
    assert result2 == expected_result_weird


@pytest.fixture
def rle1():

    return Rle([3, 9, 8, 1, 2], [5, 4, 8, 1, 3])



@pytest.fixture
def rle2():

    return Rle([8, 8, 7], [4.2, 8.0, 7.3])

@pytest.fixture
def expected_result():

    runs = [int(i) for i in "3    5    4    4    4    1    2".split()]
    values = [float(f) for f in "21 16.8   32   64 58.4  7.3 21.9".split()]

    return Rle(runs, values)

def test_mult(rle1, rle2, expected_result):

    res = rle1 * rle2

    assert res == expected_result
