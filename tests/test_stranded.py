import pyranges as pr
import numpy as np

np.random.seed(0)


def test_stranded():
    cpg = pr.data.cpg()
    exons = pr.data.exons()

    j = cpg.join(exons)

    assert j.stranded

    j.Strand = "."

    assert not j.stranded

    j.Strand = np.random.choice("+ -".split(), size=len(j))

    assert j.stranded

    for _, df in j:
        assert len(df.Strand.drop_duplicates()) == 1


def test_unstrand():

    exons = pr.data.exons()

    cpg = pr.data.cpg()
    x = exons.unstrand()
    for _, df in x:
        print(len(df.index), len(set(df.index)))
        assert len(df.index) == len(set(df.index))

    # x = pr.PyRanges(x.df.reset_index(drop=True))
    cpg.join(x)
