import pandas as pd
from sorted_nearest import makewindows


def _windows(df, kwargs):

    window_size = kwargs["window_size"]
    tile = kwargs.get("tile")

    chromosome, strand = kwargs["chromosome"], kwargs.get("strand", None)

    idxs, starts, ends = makewindows(df.index.values, df.Start.values,
                                     df.End.values, window_size, tile)

    if kwargs.get("keep_metadata"):
        df = df.reindex(idxs)
        df.loc[:, "Start"] = starts
        df.loc[:, "End"] = ends
    else:
        df = pd.DataFrame({
            "Chromosome":
            pd.Series(chromosome, dtype="category", index=idxs),
            "Start":
            starts,
            "End":
            ends
        })

        if strand:
            df.insert(3, "Strand", strand)

    return df
