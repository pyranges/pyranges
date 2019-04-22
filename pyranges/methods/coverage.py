from pyranges.multithreaded import pyrange_apply_single
from pyranges.pyranges import fill_kwargs


def _coverage(ranges, value_col=None, strand=True, rpm=False, **kwargs):

    try:
        from pyrle.methods import coverage
        from pyrle import PyRles
    except ImportError:
        raise Exception(
            "Using the coverage method requires that pyrle is installed.")

    keep = [value_col if not value_col is None else "Score"]
    kwargs = {
        "value_col": value_col,
        "sparse": {
            "self": False
        }
    }  # already sparse

    result = pyrange_apply_single(coverage, ranges, strand, kwargs)

    if rpm:
        multiplier = 1e6 / len(ranges)
        result = {k: v * multiplier for k, v in result.items()}

    return PyRles(result)
