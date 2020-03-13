def _sort(df, **kwargs):

    if "by" in kwargs:
        by = kwargs["by"]
        # print("sorting by", by)
        # print("by " * 100)
        # print(by)
        return df.sort_values(by)
    else:
        # print("else " * 100)
        # print(df)
        df = sort_one_by_one(df, "Start", "End")
        # df = df.sort_values("Start End".split())
        return df


def sort_one_by_one(d, col1, col2):
    """
    Equivalent to pd.sort_values(by=[col1, col2]), but faster.
    """

    d = d.sort_values(by=[col2])
    return d.sort_values(by=[col1], kind='mergesort')
