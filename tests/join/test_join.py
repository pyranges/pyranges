import pyranges as pr


# with slack
def test_join_with_slack():
    gr1 = pr.PyRanges(chromosomes="chr1", starts=[0], ends=[10], strands="+")
    gr2 = pr.PyRanges(chromosomes="chr1", starts=[15], ends=[20], strands="+")

    result = gr1.join(gr2, slack=10)
    df = result.df
    print(df)
    assert not df.empty
