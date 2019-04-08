import pyranges.raymock as ray

import numpy as np
import pandas as pd
from ncls import NCLS


def _coverage(ranges, value_col=None, strand=True, rpm=False):

    try:
        from pyrle.methods import coverage
        from pyrle import PyRles
    except ImportError:
        raise Exception(
            "Using the coverage method requires that pyrle is installed.")

    kwargs = {"value_col": value_col}
    # from pydbg import dbg

    result = pyrange_apply_single(coverage, ranges, strand, kwargs)

    if rpm:
        multiplier = 1e6 / len(ranges)
        result = {k: v * multiplier for k, v in result.items()}

    return PyRles(result)
