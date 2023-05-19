from io import StringIO

import numpy as np
import pandas as pd
import pytest

import pyranges as pr


def string_to_df(c):
    return pd.read_csv(
        StringIO(c),
        sep=r"\s+",
        dtype={
            "Chromosome": "category",
            "Start": np.int32,
            "End": np.int32,
            "Start_b": np.int32,
            "End_b": np.int32,
        },
    )


@pytest.fixture
def gr():
    return pr.PyRanges(chromosomes="chr1", starts=[1, 15, 200], ends=[10, 20, 2000])


@pytest.fixture
def gr2():
    return pr.PyRanges(chromosomes="chr1", starts=[11, 11, 20, 20, 50], ends=[16, 20, 21, 22, 100])


@pytest.fixture
def result_k_nearest_different():
    c = """Chromosome  Start   End  Start_b  End_b  Distance
chr1      1    10       11     16         2
chr1      1    10       11     20         2
chr1      1    10       20     21        11
chr1      1    10       20     22        11
chr1      1    10       50    100        41
chr1     15    20       11     20         0
chr1     15    20       11     16         0
chr1     15    20       20     21         1
chr1     15    20       20     22         1
chr1    200  2000       50    100      -101"""

    df = string_to_df(c)

    return df


@pytest.fixture
def result_k_nearest_different2():
    c = """Chromosome  Start  End  Start_b  End_b  Distance
0        chr1     11   16       15     20         0
1        chr1     11   20       15     20         0
2        chr1     11   20        1     10        -2
3        chr1     20   21       15     20         1
4        chr1     20   21        1     10       -11
5        chr1     20   21      200   2000       180
6        chr1     20   22       15     20         1
7        chr1     20   22        1     10       -11
8        chr1     20   22      200   2000       179
9        chr1     50  100       15     20       -31
10       chr1     50  100        1     10       -41"""

    df = string_to_df(c)

    return df


def test_k_nearest_different(result_k_nearest_different, result_k_nearest_different2, gr, gr2):
    k = [3, 2, 1]
    r = gr.k_nearest(gr2, ties="different", k=k).df

    print(r)
    print(result_k_nearest_different)

    pd.testing.assert_frame_equal(r, result_k_nearest_different)

    ks = [1, 2, 3, 4, 2]
    r2 = gr2.k_nearest(gr, ties="different", k=ks).df

    # gr2.k = ks
    print(gr2)
    print(gr)

    print("actual")
    print(r2)
    print("expected")
    print(result_k_nearest_different2)

    pd.testing.assert_frame_equal(r2, result_k_nearest_different2)


@pytest.fixture
def result_k_nearest_first():
    c = """Chromosome  Start   End  Start_b  End_b  Distance
0       chr1      1    10       11     16         2
1       chr1      1    10       20     21        11
2       chr1      1    10       50    100        41
3       chr1     15    20       11     20         0
4       chr1     15    20       20     21         1
5       chr1    200  2000       50    100      -101"""

    df = string_to_df(c)

    return df


@pytest.fixture
def result_k_nearest_first2():
    c = """Chromosome  Start  End  Start_b  End_b  Distance
0        chr1     11   16       15     20         0
1        chr1     11   20       15     20         0
2        chr1     11   20        1     10        -2
3        chr1     20   21       15     20         1
4        chr1     20   21        1     10       -11
5        chr1     20   21      200   2000       180
6        chr1     20   22       15     20         1
7        chr1     20   22        1     10       -11
8        chr1     20   22      200   2000       179
9        chr1     50  100       15     20       -31
10       chr1     50  100        1     10       -41
11       chr1     50  100      200   2000       101"""

    df = string_to_df(c)

    return df


def test_k_nearest_first(result_k_nearest_first, result_k_nearest_first2, gr, gr2):  # result_k_nearest_first,
    print(gr)
    print(gr2)
    k = [3, 2, 1]
    r = gr.k_nearest(gr2, ties="first", k=k).df

    # print(r)
    # print(result_k_nearest_first)
    # assert 0

    pd.testing.assert_frame_equal(r, result_k_nearest_first)

    r2 = gr2.k_nearest(gr, ties="first", k=list(range(1, 5 + 1))).df

    # print(r2)
    # print(result_k_nearest_first2)
    # raise

    pd.testing.assert_frame_equal(r2, result_k_nearest_first2)


@pytest.fixture
def result_k_nearest_last():
    c = """Chromosome  Start   End  Start_b  End_b  Distance
0       chr1      1    10       11     20         2
1       chr1      1    10       20     22        11
2       chr1      1    10       50    100        41
3       chr1     15    20       11     20         0
4       chr1     15    20       20     22         1
5       chr1    200  2000       50    100      -101"""

    df = string_to_df(c)

    return df


@pytest.fixture
def result_k_nearest_last2():
    c = """Chromosome  Start  End  Start_b  End_b  Distance
0        chr1     11   16       15     20         0
1        chr1     11   20       15     20         0
2        chr1     11   20        1     10        -2
3        chr1     20   21       15     20         1
4        chr1     20   21        1     10       -11
5        chr1     20   21      200   2000       180
6        chr1     20   22       15     20         1
7        chr1     20   22        1     10       -11
8        chr1     20   22      200   2000       179
9        chr1     50  100       15     20       -31
10       chr1     50  100        1     10       -41
11       chr1     50  100      200   2000       101"""

    df = string_to_df(c)

    return df


def test_k_nearest_last(result_k_nearest_last, result_k_nearest_last2, gr, gr2):  # result_k_nearest_last,
    # print(gr)
    # print(gr2)
    k = [3, 2, 1]
    r = gr.k_nearest(gr2, ties="last", k=k).df

    print(gr)
    print(gr2)
    print(r)
    print(result_k_nearest_last)

    pd.testing.assert_frame_equal(r, result_k_nearest_last)

    r2 = gr2.k_nearest(gr, ties="last", k=[1, 2, 3, 4, 5]).df

    # print("results")
    print(gr2)
    print(gr)
    print(r2)
    print(result_k_nearest_last2)

    pd.testing.assert_frame_equal(r2, result_k_nearest_last2)
    # assert 0
