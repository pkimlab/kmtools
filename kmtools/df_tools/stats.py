import numpy as np
import scipy as sp


def pearsonr_rows(df_values):
    """Calculate correlation between every row in :class:`numpy.Array` `df_values`.
    """
    nrows = df_values.shape[0]
    results = np.zeros((nrows, 2))
    for i in range(nrows):
        for j in range(i + 1, nrows):
            results[i, :] = sp.stats.pearsonr(df_values[i, :], df_values[j, :])
    return results
