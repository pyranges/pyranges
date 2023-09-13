import pandas as pd
from sorted_nearest import nearest_nonoverlapping  # type: ignore
from sorted_nearest import nearest_next_nonoverlapping, nearest_previous_nonoverlapping

from pyranges.methods.join import _both_dfs
from pyranges.methods.sort import sort_one_by_one


def _insert_distance(ocdf, dist, suffix):
    if "Distance" not in ocdf:
        distance_column_name = "Distance"
    elif "Distance" + suffix not in ocdf:
        distance_column_name = "Distance" + suffix
    else:
        i = 1
        while "Distance" + str(i) in ocdf:
            i += 1
        distance_column_name = "Distance" + str(i)

    ocdf.insert(
        ocdf.shape[1],
        distance_column_name,
        pd.Series(dist, index=ocdf.index).fillna(-1).astype(int),
    )

    return ocdf


def _overlapping_for_nearest(scdf, ocdf, suffix):
    nearest_df = pd.DataFrame(columns="Chromosome Start End Strand".split())

    scdf2, ocdf2 = _both_dfs(scdf, ocdf, how="first")

    if not ocdf2.empty:
        original_idx = scdf.index

        idxs = scdf2.index
        original_idx = scdf.index.copy(deep=True)
        missing_idxs = ~original_idx.isin(idxs)
        missing_overlap = scdf.index[missing_idxs]

        df_to_find_nearest_in = scdf.reindex(missing_overlap)

        odf = ocdf.reindex(ocdf2.index)
        odf.index = idxs
        sdf = scdf.reindex(idxs)

        nearest_df = sdf.join(odf, rsuffix=suffix)
        nearest_df = _insert_distance(nearest_df, 0, suffix)
    else:
        df_to_find_nearest_in = scdf

    return nearest_df, df_to_find_nearest_in


def _next_nonoverlapping(left_ends, right_starts, right_indexes):
    left_ends = left_ends.sort_values()
    right_starts = right_starts.sort_values()
    r_idx, dist = nearest_next_nonoverlapping(left_ends.values - 1, right_starts.values, right_indexes)
    r_idx = pd.Series(r_idx, index=left_ends.index).sort_index().values
    dist = pd.Series(dist, index=left_ends.index).sort_index().values

    return r_idx, dist


def _previous_nonoverlapping(left_starts, right_ends):
    left_starts = left_starts.sort_values()
    right_ends = right_ends.sort_values()
    r_idx, dist = nearest_previous_nonoverlapping(left_starts.values, right_ends.values - 1, right_ends.index.values)

    r_idx = pd.Series(r_idx, index=left_starts.index).sort_index().values
    dist = pd.Series(dist, index=left_starts.index).sort_index().values

    return r_idx, dist


def _nearest(scdf, ocdf, **kwargs):
    if scdf.empty or ocdf.empty:
        return None

    overlap = kwargs["overlap"]
    how = kwargs["how"]
    suffix = kwargs["suffix"]

    if how == "upstream":
        strand = scdf.Strand.iloc[0]
        how = {"+": "previous", "-": "next"}[strand]
    elif how == "downstream":
        strand = scdf.Strand.iloc[0]
        how = {"+": "next", "-": "previous"}[strand]

    ocdf = ocdf.reset_index(drop=True)

    if overlap:
        nearest_df, df_to_find_nearest_in = _overlapping_for_nearest(scdf, ocdf, suffix)
    else:
        df_to_find_nearest_in = scdf

    if not df_to_find_nearest_in.empty:
        df_to_find_nearest_in = sort_one_by_one(df_to_find_nearest_in, "Start", "End")
        ocdf = sort_one_by_one(ocdf, "Start", "End")
        df_to_find_nearest_in.index = pd.Index(range(len(df_to_find_nearest_in)))

        if how == "next":
            r_idx, dist = _next_nonoverlapping(df_to_find_nearest_in.End, ocdf.Start, ocdf.index.values)
        elif how == "previous":
            r_idx, dist = _previous_nonoverlapping(df_to_find_nearest_in.Start, ocdf.End)
        else:
            previous_r_idx, previous_dist = _previous_nonoverlapping(df_to_find_nearest_in.Start, ocdf.End)

            next_r_idx, next_dist = _next_nonoverlapping(df_to_find_nearest_in.End, ocdf.Start, ocdf.index.values)

            r_idx, dist = nearest_nonoverlapping(previous_r_idx, previous_dist, next_r_idx, next_dist)

        ocdf = ocdf.reindex(r_idx)

        ocdf.index = df_to_find_nearest_in.index

        ocdf = _insert_distance(ocdf, dist, suffix)

        r_idx = pd.Series(r_idx, index=ocdf.index)
        df_to_find_nearest_in = df_to_find_nearest_in.drop(r_idx.loc[r_idx == -1].index)

        df = df_to_find_nearest_in.join(ocdf, rsuffix=suffix)

    if overlap and "df" in locals() and not df.empty and not nearest_df.empty:
        df = pd.concat([nearest_df, df], sort=False)
    elif overlap and not nearest_df.empty:
        df = nearest_df

    df = df.drop("Chromosome" + suffix, axis=1)
    return df
