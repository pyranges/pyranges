import pyranges as pr


def test_unstrand():
    gr = pr.data.chipseq()
    u = gr.unstrand()
    assert all([not isinstance(k, tuple) for k in u.dfs])
    assert not u.stranded
