# import pyranges as pr
# import pandas as pd

# import pytest


# def test_new_position():

#     gr = pr.data.chipseq()
#     gr2 = pr.data.chipseq_background()

#     ju = gr.join(gr2, new_pos="union", suffixes="_A _B".split()).df

#     j = gr.join(gr2, suffix="_B")
#     ju2 = j.new_position("union", suffixes=(" ", "_B")).df
#     ju2.columns = ju.columns  # check_names did not work

#     pd.testing.assert_frame_equal(ju, ju2, check_names=False)
