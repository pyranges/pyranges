

# if start - ostart < 0:
#     rowdicts.append({"Chromosome": chromosome, "Start":
#                      start, "End": ostart, "Index": idx})
# elif end - oend > 0:
#     rowdicts.append({"Chromosome": chromosome, "Start":
#                         end, "End": oend - 1, "Index": idx})

# if not hits:
#     rowdicts.append({"Chromosome": chromosome, "Start": start,
#                         "End": end, "Index": idx})

# def inverse_intersection(self, other, strandedness=None):

#     sidx, oidx = both_indexes(self, other, strandedness)

#     sdf = self.df.loc[sidx]
#     odf = other.df.loc[oidx, ["Start", "End"]]

#     # max(start, ostart), "End": min(end, oend),
#     sdf.loc["Start"] = np.where(sdf.Start > odf.Start, sdf.Start, odf.Start)
#     sdf.loc["End"] = np.where(sdf.End < odf.End, sdf.End, odf.End)

#     return sdf
