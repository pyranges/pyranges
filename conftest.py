import pandas as pd
import pytest
from _pytest.doctest import DoctestItem

# The doctests show pandas 3 output, where string columns have the "str" dtype.
# pandas 2.3 can emulate this with the future.infer_string option; older pandas cannot.
_PANDAS_VERSION = tuple(int(v) for v in pd.__version__.split(".")[:2])
_EMULATE_STR_DTYPE = (2, 3) <= _PANDAS_VERSION < (3, 0)


def pytest_collection_modifyitems(config, items):
    if _PANDAS_VERSION >= (2, 3):
        return
    skip = pytest.mark.skip(reason="doctest output requires pandas >= 2.3")
    for item in items:
        if isinstance(item, DoctestItem):
            item.add_marker(skip)


def pytest_runtest_setup(item):
    if _EMULATE_STR_DTYPE and isinstance(item, DoctestItem):
        pd.set_option("future.infer_string", True)


def pytest_runtest_teardown(item, nextitem):
    if _EMULATE_STR_DTYPE and isinstance(item, DoctestItem):
        pd.set_option("future.infer_string", False)
