from __future__ import print_function

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

import pyranges.genomicfeatures.genomicfeatures as gf

random = gf.random

from pyranges.methods.itergrs import itergrs
iter = itergrs

from pyranges import statistics
stats = statistics
