from pyranges.multithreaded import pyrange_apply_single


def _to_rle(ranges, value_col=None, strand=True, rpm=False, **kwargs):
    try:
        from pyrle import PyRles  # type: ignore
        from pyrle.methods import coverage  # type: ignore
    except ImportError:
        raise Exception("Using the coverage method requires that pyrle is installed.")

    _kwargs = {
        "strand": strand,
        "value_col": value_col,
        "sparse": {"self": False},
    }  # already sparse
    kwargs.update(_kwargs)

    result = pyrange_apply_single(coverage, ranges, **kwargs)

    if rpm:
        multiplier = 1e6 / len(ranges)
        result = {k: v * multiplier for k, v in result.items()}

    return PyRles(result)
