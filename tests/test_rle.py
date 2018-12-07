import pytest

from hypothesis import given, settings, reproduce_failure, unlimited, HealthCheck, seed

from tests.hypothesis_helper import runlengths, dfs_min, runlengths_same_length_integers

from itertools import product
import tempfile
import subprocess
from io import StringIO

from pyrle import Rle

import pandas as pd
import numpy as np

# using assert df equal, because we want to consider output from bedtools and
# pyranges equal even if they have different sort order
from tests.helpers import assert_df_equal

import numpy as np

from os import environ

if environ.get("TRAVIS"):
    max_examples = 100
    deadline = None
else:
    max_examples = 100
    deadline = None


rle_operations = "+ - / *".split()

rle_operation_cmd = "Rscript --vanilla tests/compute_Rle.R {} {} '{}' {}"

@pytest.mark.r
@given(runlengths=runlengths, runlengths2=runlengths)
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@pytest.mark.parametrize("operation", rle_operations)
def test_rle(runlengths, runlengths2, operation):

    # Only compared against bioc with integers because float equality is hard,
    # for both libraries, sometimes end up with slightly different runlengths
    # when consecutive values are almost equal

    pyop = {"+": "__add__", "-": "__sub__", "*": "__mul__", "/": "__truediv__"}[operation]

    print("runlengths", runlengths)
    print("runlengths2", runlengths2)

    r = Rle(runlengths.Runs, runlengths.Values)
    r2 = Rle(runlengths2.Runs, runlengths2.Values)

    print("r\n", r)
    print("r2\n", r2)

    m = getattr(r, pyop)
    result_pyranges = m(r2)

    print("pyranges result\n", result_pyranges)

    result_df = None
    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.txt".format(temp_dir)
        f2 = "{}/f2.txt".format(temp_dir)
        outfile = "{}/result.txt".format(temp_dir)
        runlengths.to_csv(f1, sep="\t", index=False)
        runlengths2.to_csv(f2, sep="\t", index=False)

        cmd = rle_operation_cmd.format(f1, f2, operation, outfile) + " 2>/dev/null"

        subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()

        result = pd.read_table(outfile)
        s4vectors_result = Rle(result.Runs, result.Values)

    print("pyranges result\n", result_pyranges)
    print("s4vectors result\n", s4vectors_result)

    assert np.allclose(result_pyranges.runs, s4vectors_result.runs, equal_nan=False)
    assert np.allclose(result_pyranges.values, s4vectors_result.values, equal_nan=True)


rle_commute_how = ["__add__", "__mul__"]

@pytest.mark.parametrize("how", rle_commute_how)
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_commutative_rles(gr, gr2, how):

    cv = gr.coverage(strand=True)
    cv2 = gr2.coverage(strand=True)

    method = getattr(cv, how)
    method2 = getattr(cv2, how)

    result = method(cv2)
    result2 = method2(cv)

    assert result == result2, "\n".join([str(e) for e in [cv, cv2, result, result2, "---" * 10]])

@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(df=runlengths_same_length_integers)
def test_inverse_div_mul_rles(df):

    """Testing with small integers, since small value floating points might lead to
mul then div not being equal to identity function because of float equality."""

    print(df)
    runlength = df.Runs.sum()

    cv = Rle(df.Runs.values, df.Values.values)

    newruns = np.random.permutation(df.Runs.values)
    print("newruns", newruns)
    cv2 = Rle(newruns, df.Values2.values)

    print("cv\n", cv)
    print("cv2\n", cv2)

    assert runlength == np.sum(cv.runs) and runlength == np.sum(cv2.runs)

    result = cv / cv2

    result2 = result * cv2

    print("result\n", result)
    print("result2\n", result2)

    assert np.all(np.equal(result2.runs, cv.runs))
    assert np.allclose(result2.values, cv.values)


@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
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
