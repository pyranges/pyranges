
import pytest

from pyranges.pyranges import PyRanges

import pandas as pd

from io import StringIO

def test_getsetattr(chip_10):

    chip_10.whatevz = range(len(chip_10))

    assert (chip_10.whatevz == pd.Series(range(len(chip_10)))).all()



def test_getsetattr_fails(chip_10):

    with pytest.raises(Exception):
        chip_10.whatevz = range(len(chip_10) + 1)


def test_getsetattr_with_str(chip_10):

    chip_10.whatevz = "whatevz"

    assert (chip_10.whatevz == pd.Series(["whatevz"] * len(chip_10))).all()
