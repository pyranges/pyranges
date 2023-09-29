import pandas as pd


def _sort(df, **kwargs):
    if "by" in kwargs:
        by = kwargs["by"]
        if type(kwargs["by"]) is str:
            by = [by]
        else:
            by = [f for f in by]  # making sure it's a list

        if "5" in by:
            i = by.index("5")
            byp = by.copy()
            byp[i : i + 1] = ["Start", "End"]
            byn = by.copy()
            byn[i : i + 1] = ["End", "Start"]
            byna = [False if j in (i, i + 1) else True for j, f in enumerate(byn)]

            return pd.concat(
                [
                    df[df.Strand == "+"].sort_values(byp, ascending=True),
                    df[df.Strand == "-"].sort_values(byn, ascending=byna),
                ]
            )

        else:
            by = kwargs["by"]
            return df.sort_values(by)
    else:
        df = sort_one_by_one(df, "Start", "End")
        return df


def sort_one_by_one(d, col1, col2):
    """
    Equivalent to pd.sort_values(by=[col1, col2]), but faster.
    """

    d = d.sort_values(by=[col2])
    return d.sort_values(by=[col1], kind="mergesort")
