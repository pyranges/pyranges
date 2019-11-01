import pyranges as pr

def test_mcc():

    gr1, gr2 = pr.data.chipseq(), pr.data.chipseq_background()
    g = pr.data.chromsizes()
    mcc = pr.stats.mcc([gr1, gr2], genome=g, labels=["chip", "bg"], strand=True)

    print(mcc)
