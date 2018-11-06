
from ncls import NCLS

def create_ncls(df):

    return NCLS(df.Start.values, df.End.values, df.index.values)


def find_overlaps(df, start, end):

    return create_ncls(df).find_overlap(start, end)
