from itertools import product

import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings

from tests.property_based.hypothesis_helper import deadline, dfs_no_min, max_examples

strandedness = [False, "same", "opposite"]

binary_methods = [
    "set_union",
    "set_intersect",
    "overlap",
    "nearest",
    "intersect",
    "subtract",
    "join",
]

unary_methods = [
    "merge",
    "sort",
    "cluster",
    "pc",
    "mpc",
    "spc",
    "drop_duplicate_positions",
    "drop",
]

method_chain = product(binary_methods, binary_methods)
# cannot start with an operation that makes pyrange unstranded and then try a stranded op
strandedness_chain = list(product(["same", "opposite"], strandedness)) + list(product(strandedness, [None]))


@pytest.mark.bedtools
@pytest.mark.parametrize("strandedness_chain,method_chain", product(strandedness_chain, method_chain))
@settings(
    max_examples=max_examples,
    deadline=deadline,
    print_blob=True,
    suppress_health_check=HealthCheck.all(),
)
@given(gr=dfs_no_min(), gr2=dfs_no_min(), gr3=dfs_no_min())  # pylint: disable=no-value-for-parameter
# @reproduce_failure('5.5.4', b'AXicY2RAA4xIJCoLygcAALIABg==') # test_three_in_a_row[strandedness_chain122-method_chain122]
# @reproduce_failure('5.5.4', b'AXicY2QAAUYGKGBkxM9nAAABEAAJ') # test_three_in_a_row[strandedness_chain45-method_chain45]
# @reproduce_failure('5.5.4', b'AXicY2RAA4xIJDY+AAC2AAY=') # test_three_in_a_row[strandedness_chain24-method_chain24]
def test_three_in_a_row(gr, gr2, gr3, strandedness_chain, method_chain):
    print(method_chain)

    s1, s2 = strandedness_chain
    f1, f2 = method_chain

    suffix_methods = ["nearest", "join"]

    if f1 in suffix_methods and f2 in suffix_methods:
        m1 = getattr(gr, f1)
        gr2 = m1(gr2, strandedness=s1)
        if len(gr2) > 0:
            assert gr2.Start.dtype == np.int64
            assert (gr2.Start >= 0).all() and (gr2.End >= 0).all()
        m2 = getattr(gr2, f2)
        gr3 = m2(gr3, strandedness=s2, suffix="_c")
        print(gr3)
        if len(gr3) > 0:
            assert gr3.Start.dtype == np.int64
            assert (gr3.Start >= 0).all() and (gr3.End >= 0).all()

    else:
        m1 = getattr(gr, f1)
        gr2 = m1(gr2, strandedness=s1)
        if len(gr2) > 0:
            assert gr2.Start.dtype == np.int64
            assert (gr2.Start >= 0).all() and (gr2.End >= 0).all()
        m2 = getattr(gr2, f2)
        gr3 = m2(gr3, strandedness=s2)
        print(gr3)
        if len(gr3) > 0:
            assert gr3.Start.dtype == np.int64
            assert (gr3.Start >= 0).all() and (gr3.End >= 0).all()


# @pytest.mark.bedtools
# @pytest.mark.parametrize("strandedness_chain,method_chain",
#                          product(strandedness_chain, method_chain))
# @settings(
#     max_examples=max_examples,
#     deadline=deadline,
#     suppress_health_check=HealthCheck.all())
# @given(gr=dfs_no_min(), gr2=dfs_no_min(), gr3=dfs_no_min())  # pylint: disable=no-value-for-parameter
# def test_three_in_a_row(gr, gr2, gr3, strandedness_chain, method_chain):

#     s1, s2 = strandedness_chain
#     f1, f2 = method_chain

#     # print(s1, s2)
#     # print(f1, f2)

#     m1 = getattr(gr, f1)
#     gr2 = m1(gr2, strandedness=s1)
#     m2 = getattr(gr2, f2)
#     gr3 = m2(gr3, strandedness=s2)
