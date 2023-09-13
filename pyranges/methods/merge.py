import pandas as pd
from sorted_nearest import find_clusters, merge_by  # type: ignore


def _merge(df, **kwargs):
    if df.empty:
        return None

    slack = kwargs.get("slack", 0)
    chromosome, strand = kwargs["chromosome"], kwargs.get("strand", None)

    cdf = df.sort_values("Start")

    starts, ends, number = find_clusters(cdf.Start.values, cdf.End.values, slack)

    nidx = pd.Index(range(len(starts)))
    if strand:
        cluster_df = pd.DataFrame(
            {
                "Chromosome": pd.Series(chromosome, dtype="category", index=nidx),
                "Start": starts,
                "End": ends,
                "Strand": pd.Series(strand, dtype="category", index=nidx),
            }
        )
    else:
        cluster_df = pd.DataFrame(
            {
                "Chromosome": pd.Series(chromosome, dtype="category", index=nidx),
                "Start": starts,
                "End": ends,
            }
        )

    if kwargs["count"]:
        cluster_df.insert(cluster_df.shape[1], kwargs["count_col"], number)

    return cluster_df


def _merge_by(df, **kwargs):
    if df.empty:
        return None

    slack = kwargs.get("slack", 0)
    by = kwargs["by"]
    chromosome, strand = kwargs["chromosome"], kwargs.get("strand", None)

    if isinstance(by, str):
        cdf = df.sort_values([by, "Start"])
    else:
        cdf = df.sort_values(by + ["Start"])

    if isinstance(by, str):
        new_ids = (cdf[by] != cdf[by].shift()).cumsum()
    else:
        new_ids = (cdf[by] != cdf[by].shift()).any(axis=1).cumsum()

    new_to_old = pd.DataFrame(data={"new": new_ids.drop_duplicates()})

    new_to_old = pd.concat([new_to_old, cdf.reindex(new_to_old.index)[by]], axis=1)

    cdf.insert(cdf.shape[1], "ClusterBy", new_ids)

    ids, starts, ends, number = merge_by(cdf.Start.values, cdf.End.values, cdf.ClusterBy.values, slack)

    nidx = pd.Index(range(len(starts)))

    cluster_df = pd.DataFrame(
        {
            "Chromosome": pd.Series(chromosome, dtype="category", index=nidx),
            "Start": starts,
            "End": ends,
            "by": ids,
        }
    )

    if strand:
        cluster_df.insert(
            cluster_df.shape[1],
            "Strand",
            pd.Series(strand, dtype="category", index=nidx),
        )

    cluster_df = cluster_df.merge(new_to_old, left_on="by", right_on="new")
    cluster_df = cluster_df.drop(["by", "new"], axis=1)
    cluster_df = cluster_df.rename(columns={"old": by})

    if kwargs["count"]:
        cluster_df.insert(cluster_df.shape[1], kwargs["count_col"], number)

    return cluster_df
