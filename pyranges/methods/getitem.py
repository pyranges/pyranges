import pandas as pd
import numpy as np

from pyranges.subset import (get_string, get_tuple, get_slice, get_booldict)
from pyranges.methods.drop import _drop, _keep

from pyranges import PyRanges


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
        grpby = "Chromosome" if not self.stranded else ["Chromosome", "Strand"]
        to_grpby = [self.Chromosome] if not self.stranded else [self.Chromosome, self.Strand]
        d = {k: v.iloc[:, -1] for k, v in pd.concat(to_grpby + [pd.Series(val)], axis=1).groupby(grpby)}
        dfs = get_booldict(self, d)
    else:
        raise Exception("Not valid subsetter: {}".format(str(val)))

    gr = PyRanges(dfs)
    return gr
