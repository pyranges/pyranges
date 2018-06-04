import pytest
from hypothesis import given
from hypothesis.extra.pandas import data_frames, columns, range_indexes
import hypothesis.strategies as st

from itertools import product

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
df_minsize = 25
large_dfs = data_frames(index=range_indexes(min_size=df_minsize),
                        columns=columns("Chromosome Start End Strand".split(), dtype=int),
                        rows=st.tuples(chromosomes, positions, positions, strands).map(mysort))

import pyranges as pr
import pandas as pd


df_large_or_small = st.one_of(dfs, large_dfs)

strandedness = [False, "same", "opposite"]
methods = ["set_intersection", "subtraction", "join", "intersection", "overlap"]

@settings(deadline=300)
@pytest.mark.parametrize("strandedness,method", product(strandedness, methods))
@given(df=dfs, df2=dfs)
def test_methods(df, df2, strandedness, method):

    gr = pr.PyRanges(df)
    gr2 = pr.PyRanges(df2)

    method = getattr(gr, method)
    result = method(gr2, strandedness=strandedness)

    assert isinstance(result, pr.PyRanges)

how = [False, "nonoverlapping", "previous_nonoverlapping", "next_nonoverlapping", "next", "previous"]

@settings(deadline=300)
@pytest.mark.parametrize("strandedness,how", product(strandedness, how))
@given(df=dfs, df2=dfs)
def test_nearest_methods_dont_fail(df, df2, strandedness, how):

    gr = pr.PyRanges(df)
    gr2 = pr.PyRanges(df2)

    result = gr.nearest(gr2, strandedness=strandedness, how=how)

    assert isinstance(result, pr.PyRanges)
