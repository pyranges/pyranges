import pandas as pd

def _split(df, **kwargs):

    """Blah blah.

    Example
    -------

    >>> 5
    5

    TODO: would be sped up if used numpy instead of pandas for everything?
    """

    strand = kwargs.get("strand", False)

    dtype = df.Start.dtype

    starts = df.Start
    ends = df.End
    points = pd.concat([starts, ends]).sort_values().drop_duplicates()

    _ends = points.shift(-1)

    points = points[:-1]
    _ends = _ends[:-1]
    features = pd.concat([points, _ends], axis=1).astype(dtype)
    features.columns = "Start End".split()

    features.insert(0, "Chromosome", df.Chromosome.iloc[0])
    if strand and "Strand" in df:
        features.insert(features.shape[1], "Strand", df["Strand"].iloc[0])

    return features
