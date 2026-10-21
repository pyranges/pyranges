from io import StringIO

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import pyranges as pr


@pytest.mark.parametrize("int64", [False, True])
def test_window_preserves_intervals_strands_and_metadata(int64):
    gr = pr.PyRanges(
        chromosomes=["chr1"] * 4,
        starts=[3, 10, 20, 30],
        ends=[8, 14, 21, 30],
        strands=["+"] * 4,
        int64=int64,
    )
    gr.Name = ["partial", "exact", "short", "empty"]
    original = gr.df.copy(deep=True)

    result = gr.window(2).df

    assert list(result.Start) == [3, 5, 7, 10, 12, 20, 30]
    assert list(result.End) == [5, 7, 8, 12, 14, 21, 30]
    assert list(result.Name) == ["partial"] * 3 + ["exact"] * 2 + ["short", "empty"]
    assert list(result.Strand) == ["+"] * 7
    assert result.Start.dtype == original.Start.dtype
    assert_frame_equal(gr.df, original)


@pytest.mark.parametrize("window_size", [0, -1])
def test_window_rejects_nonpositive_size(window_size):
    gr = pr.PyRanges(chromosomes=["chr1"], starts=[3], ends=[8])
    with pytest.raises(ValueError, match="window_size must be positive"):
        gr.window(window_size)


@pytest.fixture
def expected_result():
    c = """Chromosome  Start  End       Name  Score Strand
0       chr1      0   10  interval1      0      +
1       chr1      0   10  interval3      0      +
2       chr1      0   10  interval2      0      -"""

    df = pd.read_csv(StringIO(c), sep=r"\s+")

    return df


def test_windows():
    f1 = pr.data.f1()

    print(f1)

    result = f1.tile(2)

    print(result)

    df = result.df

    assert list(df.Start) == [2, 4, 8, 4, 6]
    assert list(df.End) == [4, 6, 10, 6, 8]


def test_windows2():
    c = """Chromosome  Start    End  Count
0       chr1  10200  10400      7
1       chr1  10400  10600      7
2       chr1  51400  51600      1
3       chr1  51600  51800      3
4       chr1  51800  52000      3"""

    df = pd.read_csv(StringIO(c), sep=r"\s+", nrows=5)
    # df.End -= 1
    gr = pr.PyRanges(df)

    print(gr)
    result = gr.tile(200)

    print(result)

    df = result.df

    assert list(df.Start) == [10200, 10400, 51400, 51600, 51800]
    # assert list(df.End) == [4, 6, 8, 10, 6, 8]
