import pytest

import pyranges as pr
import pandas as pd

from io import StringIO

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
