import pyranges as pr


# 3' and 5'
def test_slack():
    gr = pr.PyRanges(chromosomes="chr1", starts=[15, 300], ends=[20, 305], strands="+ -".split())
    print(gr)
    gr = gr.slack({"5": 10, "3": 5})
