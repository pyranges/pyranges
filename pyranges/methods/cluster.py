
import pandas as pd

from sorted_nearest import annotate_clusters, cluster_by


def _cluster(df, kwargs):

    if df.empty:
        return None


    slack = kwargs.get("slack", 0)

    cdf = df.sort_values("Start")

    ids = annotate_clusters(cdf.Start.values, cdf.End.values, slack)

    cdf.insert(df.shape[1], "Cluster", ids)
    return cdf

def _cluster_by(df, kwargs):

    if df.empty:
        return None

    slack = kwargs.get("slack", 0)

    by = kwargs["by"]
    cdf = df.sort_values(by)

    new_ids = (cdf[by] != cdf[by].shift()).cumsum()
    cdf.insert(cdf.shape[1], "ClusterBy", new_ids)

    cdf = cdf.sort_values(["ClusterBy", "Start"])

    ids = cluster_by(cdf.Start.values, cdf.End.values, cdf.ClusterBy.values, slack)

    cdf = cdf.drop("ClusterBy", axis=1)
    cdf.insert(cdf.shape[1], "Cluster", ids)

    return cdf
