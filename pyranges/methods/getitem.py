import numpy as np
import pandas as pd

from pyranges import PyRanges
from pyranges.methods.drop import _keep
from pyranges.subset import get_booldict, get_slice, get_string, get_tuple


def _getitem(self, val):
    if isinstance(val, list):
        dfs = _keep(self, keep=val).dfs
    elif isinstance(val, str):
        dfs = get_string(self, val)
    elif isinstance(val, tuple):
        dfs = get_tuple(self, val)
    elif isinstance(val, slice):
        dfs = get_slice(self, val)
    elif isinstance(val, dict):
        dfs = get_booldict(self, val)
    elif (isinstance(val, (pd.Series, np.ndarray))) and val.dtype == "bool":
        assert len(val) == len(self), "Boolean indexer must be same length as pyrange!"
        _length = 0
        if isinstance(val, pd.Series):
            val = val.values

        dfs = {}
        for k, df in self:
            length = len(df)
            _bool = val[_length : (length + _length)]
            dfs[k] = df[_bool]
            _length += length
    else:
        raise Exception("Not a valid subsetter: {}".format(str(val)))

    gr = PyRanges(dfs)
    return gr
