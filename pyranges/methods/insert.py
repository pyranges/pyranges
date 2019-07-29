# from pyranges.methods.join import _both_dfs

# import numpy as np
# import pandas as pd

# def _insert(scdf, ocdf, kwargs):

#     how = kwargs.get("how")
#     col = kwargs.get("columns")
#     suffix = kwargs.get("suffix")
#     overlap_only = kwargs.get("overlap_only")

#     if scdf.empty or ocdf.empty:
#         return None

#     _scdf, ocdf = _both_dfs(scdf, ocdf, how=how)

#     nix = pd.Index(range(len(_scdf)))
#     _scdf.index = nix
#     ocdf.index = nix

#     ocdf = ocdf[col]

#     df = _scdf.join(ocdf, rsuffix=suffix)

#     if not overlap_only:
#         # need to insert nan for those rows in self which had no overlap

#         scdf = scdf[~scdf.index.isin(_scdf.index)]

#         if isinstance(col, str):
#             scdf.insert(scdf.shape[1], col, np.nan)
#         else:
#             for c in col:
#                 scdf.insert(scdf.shape[1], c, np.nan)

#         df = pd.concat([df, scdf])

#     df = df.sort_values("Start End".split())

#     return df
