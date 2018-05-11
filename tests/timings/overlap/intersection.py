
from time import time
import datetime

from pyranges import GRanges

import pandas as pd



chip = "/mnt/scratch/endrebak/genomes/chip/UCSD.Aorta.H3K27me3.STL003.bed.gz"
background = "/mnt/scratch/endrebak/genomes/chip/UCSD.Aorta.Input.STL002.bed.gz"

nrows = int(1e7) # None # int(1e6)

print("starting to read")
c = pd.read_table(chip, sep="\t", usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split(), nrows=nrows)
b = pd.read_table(background, sep="\t", usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split(), nrows=nrows)

print("done reading")
start = time()
c_gr = GRanges(c)

print("first range finished")

end = time()

b_gr = GRanges(b)

print("second range finished")
end2 = time()

first = datetime.datetime.fromtimestamp(end - start).strftime('%M\t%S\n')
second = datetime.datetime.fromtimestamp(end2 - end).strftime('%M\t%S\n')


start_overlap = time()

o_gr = c_gr.intersection(b_gr)

end_overlap = time()

overlap_time = datetime.datetime.fromtimestamp(end_overlap - start_overlap).strftime('%M\t%S\n')

print(first, second, overlap_time, sep="\n")

print(o_gr)

# all reads:
# 05	26

# +--------------+----------+----------+----------+
# | Chromosome   | Start    | End      | Strand   |
# |--------------+----------+----------+----------|
# | chr1         | 10151    | 10152    | -        |
# | chr1         | 10151    | 10177    | -        |
# | chr1         | 10154    | 10177    | -        |
# | ...          | ...      | ...      | ...      |
# | chrY         | 59363278 | 59363309 | +        |
# | chrY         | 59363280 | 59363294 | +        |
# | chrY         | 59363280 | 59363309 | +        |
# +--------------+----------+----------+----------+
# GRanges object has 66589595 sequences from 25 chromosomes.
