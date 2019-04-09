from pyranges.subset import get_string, get_tuple, get_slice

from pyranges import PyRanges


def _getitem(self, val):
    if isinstance(val, str):
        df = get_string(self, val)
    elif isinstance(val, tuple):
        df = get_tuple(self, val)
    elif isinstance(val, slice):
        df = get_slice(self, val)
    else:
        raise Exception("Not valid subsetter: {}".format(str(val)))

    if not df is None:
        return PyRanges(df)
    else:
        return PyRanges({})
