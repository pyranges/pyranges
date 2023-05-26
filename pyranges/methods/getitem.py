import numpy as np
import pandas as pd

import pyranges as pr
from pyranges.methods.drop import _keep
from pyranges.subset import get_2_tuple, get_booldict, get_chromosome_strand_loc, get_slice, get_string


def _getitem(self, val):
    if isinstance(val, list):
        dfs = _keep(self, keep=val).dfs
    elif isinstance(val, str):
        dfs = get_string(self, val)
    elif isinstance(val, tuple):
        if len(val) == 2:
            dfs = get_2_tuple(self, val[0], val[1])
        elif len(val) == 3:
            dfs = get_chromosome_strand_loc(self, val[0], val[1], val[2])
        else:
            raise ValueError("Indexing tuple must be of length 2 or 3. Tuple was: {}".format(str(val)))
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

    return pr.from_dfs(dfs)
