def _spliced_subseq(scdf, **kwargs):
    if scdf.empty:
        return None

    scdf = scdf.copy()
    orig_order = scdf.index.copy()

    by = kwargs.get("by") if kwargs.get("by") else "__i__"
    if not type(by) is list:
        by = [by]

    strand = kwargs.get("strand")

    # at this point, strand is False if 1. spliced_subsequence was called with strand=False or
    #                                   2. it was called with strand=None and self is not stranded
    # or it can be '+' or '-' if:
    #  1. it was input as True to spliced_subsequence and passed  to pyrange_apply_single as True,
    #     which updates it to '-' or '+' before calling _spliced_subseq, or
    #  2. it was called with strand=None and self is stranded

    if strand:
        assert "Strand" in scdf, "Cannot have strand=True on unstranded pyranges!"

    scdf.insert(scdf.shape[1], "__length__", scdf.End - scdf.Start)
    scdf.insert(scdf.shape[1], "__i__", scdf.index)

    g = scdf.groupby(by, dropna=False)
    scdf.insert(scdf.shape[1], "__cumsum__", g.__length__.cumsum())

    start = kwargs.get("start") if kwargs.get("start") else 0
    end = kwargs.get("end") if kwargs.get("end") else scdf.__cumsum__.max()

    minstart_idx = g.__i__.first()

    if start < 0 or (end is not None and end < 0):
        # len_per_transc is total sum of exon length per transcript
        len_per_transc = scdf.loc[g.__i__.last(), by + ["__cumsum__"]].rename(columns={"__cumsum__": "__totexonlen__"})

        # exp_len_per_transc has same rows of scdf with total sum of exon length
        # had to add bits to keep the order of rows right, or merge would destroy it
        if kwargs.get("by"):
            exp_len_per_transc = (
                scdf.loc[:, by + ["__i__"]].merge(len_per_transc, on=by).set_index("__i__").loc[scdf.index]
            )
        else:
            exp_len_per_transc = scdf.loc[:, by].merge(len_per_transc, on=by).set_index("__i__").loc[scdf.index]

        if start < 0:
            start = exp_len_per_transc["__totexonlen__"] + start

        if end is not None and end < 0:
            end = exp_len_per_transc["__totexonlen__"] + end

    cs_start = g.__cumsum__.shift(1, fill_value=0)
    cs_start.loc[minstart_idx] = 0

    cs_end = scdf["__cumsum__"]

    # NOTE
    # here below, start is a scalar if originally provided > 0, or a Series if < 0
    #             end is a scalar if originally None or provided >0, or a Series if provided < 0
    if strand == "-":  # and use_strand:
        start_adjustments = start - cs_start
        adjust_start = start_adjustments > 0
        scdf.loc[adjust_start, "End"] -= start_adjustments[adjust_start].astype(scdf.End.dtype)

        end_adjustments = cs_end - end
        adjust_end = end_adjustments > 0
        scdf.loc[adjust_end, "Start"] += end_adjustments[adjust_end].astype(scdf.Start.dtype)
    else:
        start_adjustments = start - cs_start
        adjust_start = start_adjustments > 0
        scdf.loc[adjust_start, "Start"] += start_adjustments[adjust_start].astype(scdf.Start.dtype)

        end_adjustments = cs_end - end
        adjust_end = end_adjustments > 0
        scdf.loc[adjust_end, "End"] -= end_adjustments[adjust_end].astype(scdf.End.dtype)

    scdf = scdf.loc[orig_order]
    scdf = scdf[(scdf.Start < scdf.End)]

    return scdf.drop(["__i__", "__length__", "__cumsum__"], axis=1)
