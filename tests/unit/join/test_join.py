from pandas.testing import assert_series_equal

import pyranges as pr


# with slack
def test_join_with_slack():
    gr1 = pr.PyRanges(chromosomes="chr1", starts=[0], ends=[10], strands="+")
    gr2 = pr.PyRanges(chromosomes="chr1", starts=[15], ends=[20], strands="+")

    result = gr1.join(gr2, slack=10)
    df = result.df
    print(df)
    assert not df.empty


def test_join_without_reordering():
    f1 = pr.from_dict(
        {
            "Chromosome": ["chr1", "chr1", "chr1"],
            "Start": [3, 8, 5],
            "End": [6, 9, 7],
            "Name": ["interval1", "interval3", "interval2"],
        }
    )
    f2 = pr.from_dict({"Chromosome": ["chr1", "chr1"], "Start": [1, 6], "End": [2, 7], "Name": ["a", "b"]})

    lj = f1.join(f2, how="left", preserve_order=True)
    assert_series_equal(lj.Name, f1.Name)

    rj = f1.join(f2, how="right", preserve_order=True)
    assert_series_equal(rj.Name_b, f2.Name, check_names=False)

    oj = f1.join(f2, how="outer", preserve_order=True)
    assert list(oj.Name) == ["interval1", "interval3", "interval2", "-1"]
    assert list(oj.Name_b) == ["-1", "-1", "b", "a"]
