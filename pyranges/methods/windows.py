import numpy as np
from sorted_nearest import maketiles, makewindows  # type: ignore


def _windows(df, **kwargs):
    window_size = kwargs["window_size"]

    idxs, starts, ends = makewindows(df.index.values, df.Start.values, df.End.values, window_size)

    df = df.reindex(idxs)
    df.loc[:, "Start"] = starts
    df.loc[:, "End"] = ends

    return df


def _intersect_tile(df):
    overlap = np.minimum(df.End, df.__End__) - np.maximum(df.Start, df.__Start__)
    df.insert(df.shape[1], "TileOverlap", overlap)

    return df


def _tiles(df, **kwargs):
    overlap = kwargs.get("overlap")

    if overlap:
        df = df.copy()
        df.insert(df.shape[1], "__Start__", df.Start)
        df.insert(df.shape[1], "__End__", df.End)

    window_size = kwargs["tile_size"]

    idxs, starts, ends = maketiles(df.index.values, df.Start.values, df.End.values, window_size)

    df = df.reindex(idxs)
    df.loc[:, "Start"] = starts
    df.loc[:, "End"] = ends

    if overlap:
        df = _intersect_tile(df)

        df = df.drop(["__Start__", "__End__"], axis=1)

    return df
