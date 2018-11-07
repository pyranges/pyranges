
from ncls import NCLS

def create_ncls(df):

    return NCLS(df.Start.values, df.End.values, df.index.values)


def find_overlaps(df, start, end):

    n = create_ncls(df)

    idxes = []
    for r in n.find_overlap(start, end):
        idxes.append(r[2])

    return idxes
