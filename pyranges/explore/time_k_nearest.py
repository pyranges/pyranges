
# kernprof -l time_k_nearest.py && python -m line_profiler time_k_nearest.py.lprof
import pyranges as pr
from os.path import expanduser
nrows = int(2e5)
gr = pr.read_bed(expanduser("~/ucsc2.bed.gz"), nrows=nrows)
gr.knearest(gr, k=1)
