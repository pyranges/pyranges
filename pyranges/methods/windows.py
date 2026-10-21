import numpy as np
from sorted_nearest import maketiles  # type: ignore


def _windows(df, **kwargs):
    window_size = kwargs["window_size"]
    if window_size <= 0:
        raise ValueError("window_size must be positive")

    # sorted_nearest 0.0.41's makewindows uses the undefined Python 2 name
    # `long`. Split by row position so repeated labels preserve their metadata.
    starts = df.Start.to_numpy()
    ends = df.End.to_numpy()
    lengths = ends.astype(np.int64) - starts
    counts = np.maximum(1, (lengths - 1) // window_size + 1)
    rows = np.repeat(np.arange(len(df)), counts)
    offsets = np.arange(len(rows)) - np.repeat(np.cumsum(counts) - counts, counts)
    window_starts = starts[rows] + offsets * window_size

    df = df.iloc[rows].copy()
    df["Start"] = window_starts.astype(starts.dtype)
    df["End"] = np.minimum(window_starts + window_size, ends[rows]).astype(ends.dtype)

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
    df["Start"] = starts
    df["End"] = ends

    if overlap:
        df = _intersect_tile(df)

        df = df.drop(["__Start__", "__End__"], axis=1)

    return df
