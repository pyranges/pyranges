import pytest

from hypothesis import given, settings, HealthCheck
from hypothesis import reproduce_failure  # pylint: disable=unused-import

from itertools import product

import numpy as np

from tests.hypothesis_helper import dfs_min2, dfs_no_min

from os import environ

if environ.get("TRAVIS"):
    max_examples = 10
    deadline = None
else:
    max_examples = 100
    deadline = None

strandedness = [False, "same", "opposite"]

binary_methods = [
    "set_union", "set_intersect", "overlap", "nearest", "intersect",
    "subtract", "join"
]

unary_methods = [
    "merge", "sort", "cluster", "pc", "mpc", "spc", "drop_duplicate_positions",
    "drop"
]

method_chain = product(binary_methods, binary_methods)
# cannot start with an operation that makes pyrange unstranded and then try a stranded op
strandedness_chain = list(product(["same", "opposite"], strandedness)) + list(
    product(strandedness, [None]))


@pytest.mark.bedtools
@pytest.mark.parametrize("strandedness_chain,method_chain",
                         product(strandedness_chain, method_chain))
@settings(
    max_examples=max_examples,
    deadline=deadline,
    print_blob=True,
    suppress_health_check=HealthCheck.all())
@given(gr=dfs_no_min(), gr2=dfs_no_min(), gr3=dfs_no_min())  # pylint: disable=no-value-for-parameter
def test_three_in_a_row(gr, gr2, gr3, strandedness_chain, method_chain):

    s1, s2 = strandedness_chain
    f1, f2 = method_chain

    # print(s1, s2)
    # print(f1, f2)

    m1 = getattr(gr, f1)
    gr2 = m1(gr2, strandedness=s1)
    if len(gr2) > 0:
        assert gr2.Start.dtype == np.int64
        assert (gr2.Start >= 0).all() and (gr2.End >= 0).all()
    m2 = getattr(gr2, f2)
    gr3 = m2(gr3, strandedness=s2)
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
