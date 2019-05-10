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
    else:
        raise Exception("Not valid subsetter: {}".format(str(val)))

    gr = PyRanges(dfs)
    return gr
