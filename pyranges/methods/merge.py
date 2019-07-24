import pandas as pd

from sorted_nearest import find_clusters, merge_by


def _merge(df, kwargs):

    if df.empty:
        return None

    slack = kwargs.get("slack", 0)
    chromosome, strand = kwargs["chromosome"], kwargs.get("strand", None)

    cdf = df.sort_values("Start")

    starts, ends, number = find_clusters(cdf.Start.values, cdf.End.values, slack)

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
            pd.Series(strand, dtype="category", index=nidx),
            "Count": number
        })
    else:
        cluster_df = pd.DataFrame({
            "Chromosome":
            pd.Series(chromosome, dtype="category", index=nidx),
            "Start":
            starts,
            "End":
            ends,
            "Count": number
        })

    return cluster_df


def _merge_by(df, kwargs):

    if df.empty:
        return None

    slack = kwargs.get("slack", 0)
    by = kwargs["by"]
    chromosome, strand = kwargs["chromosome"], kwargs.get("strand", None)

    cdf = df.sort_values([by, "Start"])

    new_ids = (cdf[by] != cdf[by].shift()).cumsum()

    new_to_old = pd.DataFrame(data={"new": new_ids.drop_duplicates(), "old": cdf[by].drop_duplicates()})

    cdf.insert(cdf.shape[1], "ClusterBy", new_ids)

    ids, starts, ends, number = merge_by(cdf.Start.values, cdf.End.values, cdf.ClusterBy.values, slack)

    nidx = pd.Index(range(len(starts)))
    if strand:
        cluster_df = pd.DataFrame({
            "Chromosome":
            pd.Series(chromosome, dtype="category", index=nidx),
            "Start":
            starts,
            "End":
            ends,
            by:
            ids,
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
            ends,
            by:
            ids,
        })

    cluster_df = cluster_df.merge(new_to_old, left_on=by, right_on="new")
    cluster_df = cluster_df.drop([by, "new"], axis=1)
    cluster_df = cluster_df.rename(columns={"old": by})
    cluster_df.insert(cluster_df.shape[1], "Count", number)

    return cluster_df
