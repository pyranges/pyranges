import pandas as pd

import pyranges as pr

exons = pr.data.exons()
cpg = pr.data.cpg()
cpg.join_overlaps(exons.remove_strand()).subset(lambda df: df.CpG > 25)[["CpG"]].assign(lambda df: df.CpG % 10, "CpGDecile")[
    "chrX"
].slack(500)
from piedpiper import Debug

with Debug():
    cpg.join_overlaps(exons.remove_strand()).subset(lambda df: df.CpG > 25)[["CpG"]].assign(lambda df: df.CpG % 10, "CpGDecile")[
        "chrX"
    ].slack(500)
