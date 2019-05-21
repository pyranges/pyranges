import pandas as pd
from sorted_nearest import makewindows
from sorted_nearest import maketiles


def _windows(df, kwargs):

    window_size = kwargs["window_size"]

    strand = kwargs.get("strand", None)

    idxs, starts, ends = makewindows(df.index.values, df.Start.values,
                                     df.End.values, window_size)

    df = df.reindex(idxs)
    df.loc[:, "Start"] = starts
    df.loc[:, "End"] = ends

    if strand:
        df.insert(3, "Strand", strand)

    return df


def _tiles(df, kwargs):

    window_size = kwargs["tile_size"]

    strand = kwargs.get("strand", None)

    idxs, starts, ends = maketiles(df.index.values, df.Start.values,
                                   df.End.values, window_size)

    df = df.reindex(idxs)
    df.loc[:, "Start"] = starts
    df.loc[:, "End"] = ends

    if strand:
        df.insert(3, "Strand", strand)

    return df
