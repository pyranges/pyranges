def _sort(df, kwargs):

    return sort_one_by_one(df, "Start", "End")


def sort_one_by_one(d, col1, col2):
    """
    Equivalent to pd.sort_values(by=[col1, col2]), but faster.
    """

    d = d.sort_values(by=[col2])
    return d.sort_values(by=[col1], kind='mergesort')
