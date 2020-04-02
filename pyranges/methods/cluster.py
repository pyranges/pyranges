from sorted_nearest import annotate_clusters, cluster_by


def _cluster(df, **kwargs):

    if df.empty:
        return None

    slack = kwargs.get("slack", 0)
    count = kwargs.get("count", False)

    cdf = df.sort_values("Start")

    ids = annotate_clusters(cdf.Start.values, cdf.End.values, slack)

    cdf.insert(df.shape[1], "Cluster", ids)

    if count:
        _count = cdf.groupby("Cluster").Cluster.count()
        _count.name = "Count"
        cdf = cdf.merge(_count, on="Cluster")

    return cdf


def _cluster_by(df, **kwargs):

    if df.empty:
        return None

    slack = kwargs.get("slack", 0)
    count = kwargs.get("count", False)

    by = kwargs["by"]

    if isinstance(by, str):
        cdf = df.sort_values([by, "Start"])
    else:
        cdf = df.sort_values(by + ["Start"])

    if isinstance(by, str):
        new_ids = (cdf[by] != cdf[by].shift()).cumsum()
    else:
        new_ids = (cdf[by] != cdf[by].shift()).any(axis=1).cumsum()

    cdf.insert(cdf.shape[1], "ClusterBy", new_ids)

    cdf = cdf.sort_values(["ClusterBy", "Start"])

    ids = cluster_by(cdf.Start.values, cdf.End.values, cdf.ClusterBy.values,
                     slack)

    cdf = cdf.drop("ClusterBy", axis=1)
    cdf.insert(cdf.shape[1], "Cluster", ids)

    if count:
        _count = cdf.groupby("Cluster").Cluster.count()
        _count.name = "Count"
        cdf = cdf.merge(_count, on="Cluster")

    return cdf
