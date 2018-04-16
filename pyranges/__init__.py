import pandas as pd

from pyranges.pyranges import GRanges

def load_dataset(basename):

    df = pd.read_table("example_data/{}.bed".format(basename), header=None,
                       names="Chromosome Start End Name Score Strand".split())

    return GRanges(df)
