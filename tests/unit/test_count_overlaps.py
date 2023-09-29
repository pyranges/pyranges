import numpy as np

import pyranges as pr
from tests.helpers import assert_df_equal

a = """Chromosome Start End Strand
chr1    6    12  +
chr1    10    20 +
chr1    22    27 -
chr1    24    30 -"""
b = """Chromosome Start End Strand
chr1    12    32 +
chr1    14    30 +"""
c = """Chromosome Start End Strand
chr1    8    15 +
chr1    713800    714800 -
chr1    32    34 -"""

grs = {n: pr.from_string(s) for n, s in zip(["a", "b", "c"], [a, b, c])}
unstranded_grs = {n: gr.unstrand() for n, gr in grs.items()}

features = pr.PyRanges(
    chromosomes=["chr1"] * 4,
    starts=[0, 10, 20, 30],
    ends=[10, 20, 30, 40],
    strands=["+", "+", "+", "-"],
)
unstranded_features = features.unstrand()


def test_strand_vs_strand_same():
    expected_result = pr.from_string(
        """Chromosome Start End Strand a b c
chr1  0 10  + 1 0 1
chr1 10 20  + 2 2 1
chr1 20 30  + 0 2 0
chr1 30 40  - 0 0 1"""
    )

    res = pr.count_overlaps(grs, features, strandedness="same")
    res = res.apply(lambda df: df.astype({"a": np.int64, "b": np.int64, "c": np.int64}))

    res.print(merge_position=True)

    assert_df_equal(res.df, expected_result.df)


# def test_strand_vs_strand_opposite():

#     expected_result = pr.from_string("""Chromosome Start End Strand a b c
# chr1  0 10  + 1 0 1
# chr1 10 20  + 1 2 1
# chr1 20 30  + 0 2 0
# chr1 30 40  - 0 0 1""")

#     res = pr.count_overlaps(grs, features, strandedness="opposite")

#     print("features")
#     print(features)

#     for name, gr in grs.items():
#         print(name)
#         print(gr)

#     res.print(merge_position=True)

#     assert_df_equal(res.df, expected_result.df)
