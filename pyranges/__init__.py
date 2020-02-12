from __future__ import print_function

try:
    import mkl
    mkl.set_num_threads(1)
except ImportError:
    pass

import pandas as pd
import numpy as np

import pkg_resources

from pyranges.pyranges import PyRanges
from pyranges.readers import read_gtf, read_bam, read_bed, read_gff3
from pyranges import data
from pyranges.methods.concat import concat

from pyrle import PyRles, Rle

from pyranges.version import __version__

get_example_path = data.get_example_path

read_gff = read_gtf

def from_dict(d):
    return PyRanges(pd.DataFrame(d))

import pyranges.genomicfeatures.genomicfeatures as gf

random = gf.random

from pyranges.methods.itergrs import itergrs
iter = itergrs

from pyranges.methods.multioverlap import count_overlaps#, interval_split

from pyranges import statistics
stats = statistics
