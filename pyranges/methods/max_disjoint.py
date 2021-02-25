#!/usr/bin/env python

from sorted_nearest import max_disjoint


def _max_disjoint(df, **kwargs):

    if df.empty:
        return None

    slack = kwargs.get("slack", 0)

    cdf = df.sort_values("End")

    idx = max_disjoint(cdf.index.values, cdf.Start.values, cdf.End.values, slack)

    return cdf.reindex(idx)
