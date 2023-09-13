def _drop_duplicate_positions(df, **kwargs):
    strand = kwargs.get("strand")
    keep = kwargs.get("keep")

    columns = ["Start", "End"]
    if strand:
        columns.append("Strand")

    return df.drop_duplicates(columns, keep=keep)
