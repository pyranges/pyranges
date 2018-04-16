import pandas as pd

import pkg_resources

from pyranges.pyranges import GRanges

def load_dataset(basename):

    full_path = pkg_resources.resource_filename("pyranges", "example_data/{}.bed".format(basename))

    df = pd.read_table(full_path, header=None,
                       names="Chromosome Start End Name Score Strand".split())

    return GRanges(df)


def list_datasets():

    datasets = [f.replace(".bed", "") for f in pkg_resources.resource_listdir("pyranges", "example_data")]

    print(datasets)
