import pandas as pd
import numpy as np

from ncls import NCLS


def add_rows_per_group(df):
    last_rows = df.groupby("__ix__").last().reset_index()
    last_rows.loc[:, "__last__"] = True
    df = pd.concat([df, last_rows], ignore_index=True)
    df = df.sort_values("__ix__", ascending=True)
    return df

# def _subtraction(scdf, **kwargs):
#     if scdf.empty:
#         return scdf

#     falses = np.zeros(len(scdf), dtype=bool)
#     scdf.insert(scdf.shape[1], "__first__", falses)
#     scdf.insert(scdf.shape[1], "__last__", falses)

#     scdf = add_rows_per_group(scdf)

#     scdf.insert(scdf.shape[1], "NewStart", scdf.End__deleteme__.shift(fill_value=-1))
#     scdf.insert(scdf.shape[1], "NewEnd", scdf.Start__deleteme__)
#     scdf.insert(scdf.shape[1], "__ix2__", np.arange(len(scdf)))

#     first_rows = scdf.groupby(scdf.__ix__, as_index=False).first()
#     scdf.loc[scdf.__ix2__.isin(first_rows.__ix2__), "__first__"] = True

#     scdf.loc[:, "NewStart"] = np.where(scdf.__first__, scdf.Start, scdf.NewStart)

#     scdf.loc[scdf.__first__ & ~(scdf.Start__deleteme__ >= scdf.Start), "NewStart"] = -1
#     scdf.loc[:, "NewEnd"] = np.where(scdf.__last__, scdf.End, scdf.NewEnd)
#     scdf.loc[:, "NewStart"] = np.where(scdf.__last__, scdf.End__deleteme__, scdf.NewStart)
#     scdf.loc[scdf.__last__ & ~(scdf.End__deleteme__ <= scdf.End), ["NewEnd", "NewStart"]] = -1


#     scdf = scdf[~((scdf.NewStart == -1) | (scdf.NewEnd == -1))]
#     scdf = scdf.drop(["Start", "End"], axis=1)
#     scdf.rename(columns={"NewStart": "Start", "NewEnd": "End"}, inplace=True)

#     remove_mask = scdf.Start >= scdf.End

#     return scdf[~remove_mask]

def _subtraction(scdf, ocdf, **kwargs):

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
        scdf.Start.values, scdf.End.values, scdf.index.values, scdf.__num__.values)

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
