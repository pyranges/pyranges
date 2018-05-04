from pyranges import GRanges


if __name__ == "__main__":

    "kernprof -l pyrle/rledict.py && python -m line_profiler coverage.py.lprof"

    from time import time
    import datetime

    import pandas as pd

    test_file = "/mnt/scratch/endrebak/genomes/chip/UCSD.Aorta.Input.STL002.bed.gz"

    df = pd.read_table(test_file, sep="\t", usecols=[0, 1, 2, 5], header=None,
                       names="Chromosome Start End Strand".split(), nrows=None)


    print("Done reading")
    start = time()

    result = GRanges(df)

    end = time()
    total = end - start

    total_dt = datetime.datetime.fromtimestamp(total)

    minutes_seconds = total_dt.strftime('%M\t%S\n')

    print(result)
    print(minutes_seconds)
