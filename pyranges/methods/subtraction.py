import pandas as pd
from ncls import NCLS


def _subtraction(scdf, ocdf, kwargs):

    if ocdf.empty or scdf.empty:
        return scdf

    strandedness = kwargs["strandedness"]
    strand = True if strandedness else False

    chromosome = scdf.Chromosome.head(1).iloc[0]
    kwargs["chromosome"] = chromosome

    if "Strand" in ocdf and strand:
        strand = scdf.Strand.head(1).iloc[0]
        kwargs["strand"] = strand

    o = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    idx_self, new_starts, new_ends = o.set_difference_helper(
        scdf.Start.values, scdf.End.values, scdf.index.values)

    missing_idx = pd.Index(scdf.index).difference(idx_self)

    idx_to_drop = new_starts != -1

    new_starts = new_starts[idx_to_drop]
    new_ends = new_ends[idx_to_drop]

    idx_self = idx_self[idx_to_drop]
    new_starts = pd.Series(new_starts, index=idx_self)
    new_ends = pd.Series(new_ends, index=idx_self)

    scdf = scdf.reindex(missing_idx.union(idx_self)).sort_index()
    new_starts = new_starts.sort_index()
    new_ends = new_ends.sort_index()

    if len(idx_self):
        scdf.loc[scdf.index.isin(idx_self), "Start"] = new_starts.values
        scdf.loc[scdf.index.isin(idx_self), "End"] = new_ends.values

    if not scdf.empty:
        return scdf
    else:
        return None
