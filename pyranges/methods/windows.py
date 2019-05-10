import pandas as pd
from sorted_nearest import makewindows


def _windows(df, kwargs):

    window_size = kwargs["window_size"]
    tile = kwargs.get("tile")

    strand = kwargs.get("strand", None)

    idxs, starts, ends = makewindows(df.index.values, df.Start.values,
                                     df.End.values, window_size, tile)

    df = df.reindex(idxs)
    df.loc[:, "Start"] = starts
    df.loc[:, "End"] = ends

    if strand:
        df.insert(3, "Strand", strand)

    return df
