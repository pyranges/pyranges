import pyranges as pr
import pickle

def test_pickle():
    gr = pr.data.f1()
    pickle.dump(gr, open("hi", "wb+"))
    gr2 = pickle.load(open("hi", "rb"))
