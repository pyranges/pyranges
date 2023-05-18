import pandas as pd

import pyranges as pr

exons = pr.data.exons()
cpg = pr.data.cpg()
cpg.join(exons.unstrand()).subset(lambda df: df.CpG > 25)[["CpG"]].assign(lambda df: df.CpG % 10, "CpGDecile")[
    "chrX"
].slack(500)
from piedpiper import Debug

with Debug():
    cpg.join(exons.unstrand()).subset(lambda df: df.CpG > 25)[["CpG"]].assign(lambda df: df.CpG % 10, "CpGDecile")[
        "chrX"
    ].slack(500)
