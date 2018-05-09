import pandas as pd


def assert_df_equal(df1, df2):

    df1 = df1.sort_values("Chromosome Start End".split())
    df2 = df2.sort_values("Chromosome Start End".split())

    df1 = df1.reset_index(drop=True)
    df2 = df2.reset_index(drop=True)

    print("Actual")
    print(df1)
    print("Expected")
    print(df2)

    return df1.equals(df2)
