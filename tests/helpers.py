import pandas as pd


def assert_df_equal(df1, df2):

    if "Strand" in df1 and "Strand" in df2:
        sort_on = "Chromosome Start End Strand".split()
    else:
        sort_on = "Chromosome Start End".split()

    df1 = df1.sort_values(sort_on)
    df2 = df2.sort_values(sort_on)

    df1 = df1.reset_index(drop=True)
    df2 = df2.reset_index(drop=True)

    print("Actual")
    print(df1)
    print("Expected")
    print(df2)

    print("Actual dtypes")
    print(df1.dtypes)
    print("Expected dtypes")
    print(df2.dtypes)
    print("dtypes equal", df1.dtypes == df2.dtypes)

    print("Actual index")
    print(df1.index)
    print("Expected index")
    print(df2.index)
    print("index equal", df1.index == df2.index)

    return df1.equals(df2)
