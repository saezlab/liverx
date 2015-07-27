import numpy as np
from scipy.stats.stats import pearsonr, spearmanr


def pearson(x, y):
    mask = np.bitwise_and(np.isfinite(x), np.isfinite(y))
    cor, pvalue = pearsonr(x[mask], y[mask])
    return cor, pvalue, sum(mask)


def spearman(x, y):
    mask = np.bitwise_and(np.isfinite(x), np.isfinite(y))
    cor, pvalue = spearmanr(x[mask], y[mask])
    return cor, pvalue, sum(mask)
