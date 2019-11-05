import pytest 

import pyranges as pr

def test_k_nearest_first():

    gr = pr.PyRanges(chromosomes="chr1", starts=[1, 15, 200], ends=[10, 20, 2000])
    gr2 = pr.PyRanges(chromosomes="chr1", starts=[11, 11, 20, 20, 50], ends=[16, 20, 21, 22, 100])

    print(gr)
    print(gr2)
    print("----")
    # print(gr.k_nearest(gr2, ties="first", how="upstream", k=2))
    # print(gr.k_nearest(gr2, ties="first", how="downstream", k=2))
    # print(gr2.k_nearest(gr, ties="first", how="upstream", k=2))
    # print(gr2.k_nearest(gr, ties="first", how="downstream", k=2))
    k = [3, 2, 1]
    r = gr.k_nearest(gr2, ties="different", k=k)
    print(r.df)
    # print(gr2.k_nearest(gr, ties="first", k=5).df)

    # assert 0
