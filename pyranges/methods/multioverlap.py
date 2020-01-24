from io import StringIO
import pyranges as pr
import pandas as pd


def count_overlaps(grs, features=None, how=None, nb_cpu=1, strandedness=None):

    if features is None:
        features = pr.concat(grs.values()).split()

    from pyranges.methods.intersection import _count_overlaps

    hits_gr = {}
    for name, gr in grs.items():

        gr = gr.drop()

        res = features.apply_pair(gr, _count_overlaps, as_pyranges=False, nb_cpu=nb_cpu, strandedness=strandedness)

        setattr(features, name, res)

        setattr(features, name, getattr(features, name).fillna(0))

    return features

# if __name__

# a = """Chromosome Start End
# chr1    6    12
# chr1    10    20
# chr1    22    27
# chr1    24    30"""

# b = """Chromosome Start End
# chr1    12    32
# chr1    14    30"""

# c = """Chromosome Start End
# chr1    8    15
# chr1    10    14
# chr1    32    34"""

# a, b, c = [pr.PyRanges(pd.read_table(StringIO(x), sep="\s+")) for x in [a, b, c]]

# if __name__ == "__main__":

#     print(a)


