def _bounds(scdf, **kwargs):
    if scdf.empty:
        return None

    col_order = [c for c in scdf.columns]

    by = kwargs.get("group_by")
    if not type(by) is list:
        by = [by]

    agg_dict = kwargs.get("agg") if kwargs.get("agg") else {}
    agg_dict.update({"Start": "min", "End": "max", "Chromosome": "first"})
    if "Strand" in scdf.columns:
        agg_dict["Strand"] = "first"

    res = scdf.groupby(by).agg(agg_dict).reset_index()
    res = res.reindex(columns=[c for c in col_order if c in res.columns])

    return res
