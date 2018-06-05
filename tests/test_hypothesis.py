import pytest
from hypothesis import given, settings, reproduce_failure
from hypothesis.extra.pandas import data_frames, columns, range_indexes, column
from hypothesis.extra.numpy import arrays
import hypothesis.strategies as st

from itertools import product
import tempfile
import subprocess
from io import StringIO

def mysort(tp):

    if tp[1] == tp[2]:
        tp = (tp[0], tp[1], tp[2] + 1, tp[3])

    key = [-1, tp[1], tp[2], int(1e10)]

    return [x for _, x in sorted(zip(key, tp))]

chromosomes = st.sampled_from(["chr{}".format(str(e)) for e in list(range(1, 23)) + "X Y M".split()])

positions = st.integers(min_value=0, max_value=int(1e7))
strands = st.sampled_from("+ -".split())
dfs = data_frames(columns=columns("Chromosome Start End Strand".split(),
                                  dtype=int), rows=st.tuples(chromosomes, positions, positions,
                                                             strands).map(mysort))

df_minsize = 1
nonempty_dfs = data_frames(index=range_indexes(min_size=df_minsize),
                           columns=columns("Chromosome Start End Strand".split(), dtype=int),
                           rows=st.tuples(chromosomes, positions, positions, strands).map(mysort))

runlengths = data_frames(index=range_indexes(min_size=df_minsize),
                         columns=[column("Runs", st.integers(min_value=1, max_value=int(1e7))),
                                  column("Values", st.floats(min_value=0.0000001, max_value=int(1e7)))])

# runlengths_same_length_floats = data_frames(index=range_indexes(min_size=df_minsize),
#                                             columns=[column("Runs", st.integers(min_value=1, max_value=int(1e7))),
#                                                      column("Values", st.floats(min_value=0.001, max_value=int(1e7))),
#                                                      column("Values2", st.floats(min_value=0.001, max_value=int(1e7)))])


runlengths_same_length_integers = data_frames(index=range_indexes(min_size=df_minsize),
                                              columns=[column("Runs", st.integers(min_value=1, max_value=int(1e5))),
                                              column("Values", st.integers(min_value=1, max_value=int(1e5))),
                                              column("Values2", st.integers(min_value=1, max_value=int(1e5)))])


# runs = arrays(st.integers(min_value=1, max_value=int(1e5)), shape=)
# values = arrays(dtype=np.float)

import pyranges as pr
import pandas as pd
from pyrle import Rle
import numpy as np


strandedness = [False, "same", "opposite"]
methods = ["set_intersection", "subtraction", "join", "intersection", "overlap"]

@pytest.mark.parametrize("strandedness,method", product(strandedness, methods))
@settings(deadline=300)
@given(df=dfs, df2=dfs)
def test_methods(df, df2, strandedness, method):

    gr = pr.PyRanges(df)
    gr2 = pr.PyRanges(df2)

    method = getattr(gr, method)
    result = method(gr2, strandedness=strandedness)

    assert isinstance(result, pr.PyRanges)

how = [False, "nonoverlapping", "previous_nonoverlapping", "next_nonoverlapping", "next", "previous"]

@pytest.mark.parametrize("strandedness,how", product(strandedness, how))
@settings(deadline=300)
@given(df=dfs, df2=dfs)
def test_nearest_methods_dont_fail(df, df2, strandedness, how):

    gr = pr.PyRanges(df)
    gr2 = pr.PyRanges(df2)

    result = gr.nearest(gr2, strandedness=strandedness, how=how)

    assert isinstance(result, pr.PyRanges)

rle_commute_how = ["__add__", "__mul__"]

@pytest.mark.parametrize("how", rle_commute_how)
@settings(deadline=300)
@given(df=nonempty_dfs, df2=nonempty_dfs)
def test_commutative_rles(df, df2, how):

    cv = pr.PyRanges(df).coverage(stranded=True)
    cv2 = pr.PyRanges(df2).coverage(stranded=True)

    method = getattr(cv, how)
    method2 = getattr(cv2, how)

    result = method(cv2)
    result2 = method2(cv)

    assert result == result2, "\n".join([str(e) for e in [cv, cv2, result, result2, "---" * 10]])

rle_inverse_how = [["__add__", "__sub__"], ["__truediv__", "__mul__"]]

# @pytest.mark.parametrize("hows", rle_inverse_how)
@settings(deadline=300)
@given(df=runlengths_same_length_integers)
def test_inverse_div_mul_rles(df):

    """Testing with small integers, since small value floating points might lead to
mul then div not being equal to identity function because of float equality."""

    cv = Rle(df.Runs.values, df.Values.values)

    cv2 = Rle(np.random.permutation(df.Runs.values), df.Values2.values)

    result = cv / cv2

    result2 = result * cv2

    assert np.all(np.equal(result2.runs, cv.runs))
    assert np.allclose(result2.values, cv.values)

@settings(deadline=300)
@given(df=runlengths_same_length_integers)
def test_inverse_add_sub_rles(df):

    """Testing with small integers, since small value floating points might lead to
mul then div not being equal to identity function because of float equality."""

    cv = Rle(df.Runs.values, df.Values.values)

    cv2 = Rle(np.random.permutation(df.Runs.values), df.Values2.values)

    result = cv + cv2

    result2 = result - cv2

    assert np.all(np.equal(result2.runs, cv.runs))
    assert np.allclose(result2.values, cv.values)


# @settings(deadline=300)
# @given(df=nonempty_dfs, df2=nonempty_dfs)
# @reproduce_failure('3.57.0', b'AXicY2CAAUYoDQAAGgAC')
# def test_nearest_equal_to_bedtools(df, df2):

#     result_df = None
#     with tempfile.TemporaryDirectory() as temp_dir:
#         f1 = "{}/f1.bed".format(temp_dir)
#         f2 = "{}/f2.bed".format(temp_dir)
#         df.to_csv(f1, sep="\t", header=False, index=False)
#         df2.to_csv(f2, sep="\t", header=False, index=False)

#         cmd = "bedtools closest -d -a {} -b {}".format(f1, f2)
#         result = subprocess.check_output(cmd, shell=True).decode()

#         bedtools_df = pd.read_table(StringIO(result), header=None, squeeze=True, names="C S E St C2 S2 E2 St2 Distance".split())

#         distances_bedtools = list(bedtools_df.Distance.values)

#     gr = pr.PyRanges(df)
#     gr2 = pr.PyRanges(df2)

#     result = gr.nearest(gr2)

#     if not result.df.empty:
#         pyranges_distances = list(result.df.Distance.values)
#     else:
#         pyranges_distances = []

#     print("bedtools", distances_bedtools)
#     print("bedtools_df", bedtools_df)
#     print("pyranges", pyranges_distances)

#     assert distances_bedtools == pyranges_distances
