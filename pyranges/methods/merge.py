import pandas as pd

from sorted_nearest import find_clusters


def _merge(df, kwargs):

    if df.empty:
        return None

    chromosome, strand = kwargs["chromosome"], kwargs.get("strand", None)

    cdf = df.sort_values("Start")

    starts, ends = find_clusters(cdf.Start.values, cdf.End.values)

    nidx = pd.Index(range(len(starts)))
    if strand:
        cluster_df = pd.DataFrame({
            "Chromosome":
            pd.Series(chromosome, dtype="category", index=nidx),
            "Start":
            starts,
            "End":
            ends,
            "Strand":
            pd.Series(strand, dtype="category", index=nidx)
        })
    else:
        cluster_df = pd.DataFrame({
            "Chromosome":
            pd.Series(chromosome, dtype="category", index=nidx),
            "Start":
            starts,
            "End":
            ends
        })

    return cluster_df
