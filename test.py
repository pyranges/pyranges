from pyranges import PyRanges


if __name__ == "__main__":

    "kernprof -l pyrle/rledict.py && python -m line_profiler coverage.py.lprof"

    from time import time
    import datetime

    import pandas as pd

    # test_file = "/mnt/scratch/endrebak/genomes/chip/UCSD.Aorta.Input.STL002.bed.gz"
    # test_file2 = "/mnt/scratch/endrebak/genomes/chip/shuffled_Aorta_Input.gz"
    # test_file = "/local/home/endrebak/large_data/consensus_epigenome/epigenome_roadmap/H3K9me3/data/hg38/bed/Ovary_1_Input.bed"
    # test_file = "shuffled_input.bed"
    f1 = "/mnt/scratch/endrebak/pyranges_benchmark/data/download/chip_15000000.bed.gz"
    f2 = "/mnt/scratch/endrebak/pyranges_benchmark/data/download/input_15000000.bed.gz"

    nrows = int(1e6) # None # int(1e6)
    df = pd.read_table(f1, sep="\t", usecols=[0, 1, 2, 5], header=None,
                       names="Chromosome Start End Strand".split(), nrows=nrows)

    df2 = pd.read_table(f2, sep="\t", usecols=[0, 1, 2, 5], header=None,
                        names="Chromosome Start End Strand".split(), nrows=nrows)


    print("Done reading")
    start = time()

    a = PyRanges(df)
    b = PyRanges(df2)
    result = a.nearest(b)

    end = time()
    total = end - start

    total_dt = datetime.datetime.fromtimestamp(total)

    minutes_seconds = total_dt.strftime('%M\t%S\n')

    print(minutes_seconds)
    print(result)


# +--------------+----------+----------+----------+
# | Chromosome   | Start    | End      | Strand   |
# |--------------+----------+----------+----------|
# | chr1         | 10000    | 10196    | -        |
# | chr1         | 10070    | 10270    | -        |
# | chr1         | 10078    | 10278    | -        |
# | ...          | ...      | ...      | ...      |
# | chr2         | 31725669 | 31725869 | -        |
# | chr2         | 31725688 | 31725888 | -        |
# | chr2         | 31725778 | 31725978 | -        |
# +--------------+----------+----------+----------+
# PyRanges object has 10000000 sequences from 28 chromosomes.
# Wrote profile results to test.py.lprof
# Timer unit: 1e-06 s

# Total time: 39.6635 s
# File: /local/home/endrebak/code/pyranges/pyranges/pyranges.py
# Function: __init__ at line 25

# Line #      Hits         Time  Per Hit   % Time  Line Contents
# ==============================================================
#     25                                               @profile
#     26                                               def __init__(self, df):
#     27
#     28         1     598551.0 598551.0      1.5          df.Chromosome = df.Chromosome.astype("category")
#     29         1    2790564.0 2790564.0      7.0          df.Chromosome.cat.reorder_categories(natsorted(df.Chromosome.drop_duplicates()), inplace=True, ordered=True)
#     30
#     31                                                   # df = df.sort_values(["Chromosome", "Start", "End"])
#     32
#     33                                                   # so that we can be sure that an integer representation of chromosomes won't be promoted to floats
#     34         1    7480705.0 7480705.0     18.9          df.Chromosome = df.Chromosome.astype(str)
#     35
#     36         1    2396950.0 2396950.0      6.0          df = df.reset_index(drop=True)
#     37
#     38         1          4.0      4.0      0.0          self.df = df
#     39
#     40         1          2.0      2.0      0.0          self.__ncls__ = dict()
#     41
#     42         1         54.0     54.0      0.0          if "Strand" not in df:
#     43
#     44                                                       for chromosome, cdf in df.groupby("Chromosome"):
#     45
#     46                                                           indexes = cdf.index.values
#     47                                                           starts = cdf.Start.values
#     48                                                           ends = cdf.End.values
#     49
#     50                                                           self.__ncls__[chromosome, "+"] = NCLS(starts, ends, indexes)
#     51
#     52                                                   else:
#     53
#     54        37   23664466.0 639580.2     59.7              for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):
#     55
#     56        36        269.0      7.5      0.0                  indexes = cdf.index.values
#     57        36       5769.0    160.2      0.0                  starts = cdf.Start.values
#     58        36       4757.0    132.1      0.0                  ends = cdf.End.values
#     59
#     60        36    2721420.0  75595.0      6.9                  self.__ncls__[chromosome, strand] = NCLS(starts, ends, indexes)


# df = fread("zcat /mnt/scratch/endrebak/genomes/chip/shuffled_Aorta_Input.gz | cut -f 1-3,6")
# start.time <- Sys.time()
# gr = PyRanges(seqnames=df$V1, IRanges(start=df$V2, end=df$V3), strand=df$V4)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
