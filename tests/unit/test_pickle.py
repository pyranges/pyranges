import pickle

from pandas.testing import assert_frame_equal

import pyranges as pr


def test_pickle():
    gr = pr.data.f1()
    pickle.dump(gr, open("hi", "wb+"))
    gr2 = pickle.load(open("hi", "rb"))
    assert_frame_equal(gr.df, gr2.df)
