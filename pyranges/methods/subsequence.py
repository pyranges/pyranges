import numpy as np


def _subseq(scdf, **kwargs):
    if scdf.empty:
        return None

    scdf = scdf.copy()
    orig_order = scdf.index.copy()
    scdf.insert(0, "__i__", orig_order)

    by = kwargs.get("by") if kwargs.get("by") else "__i__"
    if not type(by) is list:
        by = [by]

    strand = kwargs.get("strand")
    # at this point, strand is False if:
    #   1. subsequence was called with strand=False or
    #   2. it was called with strand=None and self is not stranded
    # or it can be '+' or '-' if:
    #   1. it was input as True to subsequence and passed  to pyrange_apply_single as True,
    #      which updates it to '-' or '+' before calling _subseq, or
    #   2. it was called with strand=None and self is stranded
    if strand != "-":
        strand = "+"
    # now, unstranded or strand==None cases are treated like all intervals are on the + strand

    # creating j which holds the boundaries per group
    # j contains one row per group; columns: Start  End (+ by columns); indexed by __i__
    agg_dict = {"__i__": "first", "Start": "min", "End": "max"}
    for b in by:
        agg_dict[b] = "first"

    if kwargs.get("by"):
        j = scdf.groupby(by, dropna=False)[["Start", "End", "__i__"] + by].agg(agg_dict).set_index("__i__")
    else:
        j = scdf.groupby(by, dropna=False)[["Start", "End", "__i__"]].agg(agg_dict).set_index("__i__")
        j.insert(0, "__i__", j.index)
        j.index.name = None

    start = kwargs.get("start") if kwargs.get("start") else 0
    end = kwargs.get("end") if kwargs.get("end") else (j.End - j.Start).max()

    # below we add columns starts, ends to j
    # start and ends define the desired start and end, with one row per group (transcript).
    # they may be out of bounds of the interval, though (this is dealt with later)

    if (strand == "+" and start >= 0) or (strand == "-" and start < 0):
        j["starts"] = j.Start + abs(start)
    else:
        j["starts"] = j.End - abs(start)

    if (strand == "+" and end >= 0) or (strand == "-" and end < 0):
        j["ends"] = j.Start + abs(end)
    else:
        j["ends"] = j.End - abs(end)

    if strand == "-":
        j.rename(columns={"ends": "__min__", "starts": "__max__"}, inplace=True)
    else:
        j.rename(columns={"starts": "__min__", "ends": "__max__"}, inplace=True)

    # I'm maintaing the original row order
    scdf = scdf.merge(j[by + ["__min__", "__max__"]], on=by).set_index("__i__").loc[orig_order]

    # instead of simply using starts and ends as computed above, we're dealing here with potential out of bounds:
    r = scdf[~((scdf.Start >= scdf.__max__) | (scdf.End <= scdf.__min__))].copy()
    r.loc[:, "Start"] = np.maximum(r.Start, r.__min__)
    r.loc[:, "End"] = np.minimum(r.End, r.__max__)

    r = r.drop(["__min__", "__max__"], axis=1)

    return r
