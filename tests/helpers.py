import pandas as pd


def assert_df_equal(df1, df2):

    # df1.loc[:, "Start"] = df1.Start.astype(np.int64)
    # df2.loc[:, "Start"] = df1.Start.astype(np.int64)
    # df1.loc[:, "End"] = df1.End.astype(np.int64)
    # df2.loc[:, "End"] = df1.End.astype(np.int64)

    pd.options.mode.chained_assignment = None
    if "Strand" in df1 and "Strand" in df2:
        sort_on = "Chromosome Start End Strand".split()
        df1.Strand = df1.Strand.astype("object")
        df2.Strand = df2.Strand.astype("object")
    else:
        sort_on = "Chromosome Start End".split()

    if "Strand_b" in df1:
        sort_on += "Start_b End_b Strand_b".split()
        df1.Strand_b = df1.Strand_b.astype("object")
        df2.Strand_b = df2.Strand_b.astype("object")
    elif "Start_b" in df2:
        sort_on += "Start_b End_b".split()

    df1 = df1.sort_values(sort_on)
    df2 = df2.sort_values(sort_on)

    df1 = df1.reset_index(drop=True)
    df2 = df2.reset_index(drop=True)

    df1.Chromosome = df1.Chromosome.astype("object")
    df2.Chromosome = df2.Chromosome.astype("object")

    print("Actual")
    print(df1.to_csv(sep=" "))
    print("Expected")
    print(df2.to_csv(sep=" "))

    print("Actual dtypes")
    print(df1.dtypes)
    print("Expected dtypes")
    print(df2.dtypes)
    # print("dtypes Strand\n", "1",  df1.Strand.dtype, "2", df2.Strand.dtype)
    # print("dtypes Strand\n", df1.Strand.dtype == df2.Strand.dtype)
    # print("dtypes equal\n", df1.dtypes == df2.dtypes)

    print("Actual index")
    print(df1.index)
    print("Expected index")
    print(df2.index)
    print("index equal", df1.index == df2.index)

    pd.testing.assert_frame_equal(df1, df2)

    pd.options.mode.chained_assignment = "warn"
