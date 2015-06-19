import itertools as it
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from liverx import wd
from scipy.stats.stats import pearsonr
from pandas import DataFrame, read_csv, melt


def pearson(x, y):
    mask = np.bitwise_and(np.isfinite(x), np.isfinite(y))
    cor, pvalue = pearsonr(x[mask], y[mask]) if np.sum(mask) > 1 else (np.NaN, np.NaN)
    return cor, pvalue, sum(mask)


sns.set_style('white')

swath_quant = read_csv('%s/data/result_swath_v2.3.7_protein_quant.tab' % wd, sep='\t').replace(0.0, np.NaN)

# ---- Defined conditions
b6 = [c for c in swath_quant.columns if c.startswith('B6')]
s9 = [c for c in swath_quant.columns if c.startswith('S9')]

hypothesis, fdr_thres = 'H2', '0.01'

# ---- Import String network
network = read_csv('%s/files/string_mouse_network_filtered_800.txt' % wd, sep='\t', names=['p1', 'p2'], skiprows=1)
network = set(zip(*(network['p1'], network['p2'])))

# ---- Import subnetwork
subnetwork = read_csv('%s/files/network_enrichment/%s_%s_network.sif' % (wd, hypothesis, fdr_thres), header=None, sep='\t', names=['p1', 'i', 'p2'])
subnetwork = set(subnetwork['p1']).intersection(subnetwork['p2'])
subnetwork = {(p1, p2) for p1, p2 in network if p1 in subnetwork or p2 in subnetwork}
subnetwork = {(p1, p2) for p1, p2 in network if p1 in swath_quant.index and p2 in swath_quant.index}

# ---- Compute correlations
b6_cor = [(p1, p2, pearson(swath_quant.ix[p1, b6], swath_quant.ix[p2, b6])) for p1, p2 in subnetwork]
b6_cor = [(p1, p2, cor, pvalue, meas) for p1, p2, (cor, pvalue, meas) in b6_cor]
b6_cor = DataFrame(b6_cor, columns=['p1', 'p2', 'cor', 'pvalue', 'meas']).dropna()
b6_cor = b6_cor[b6_cor['meas'] > 2]
b6_cor['adj.pvalue'] = multipletests(b6_cor['pvalue'], method='fdr_bh')[1]
b6_cor = b6_cor.sort('adj.pvalue')
b6_cor['strain'] = 'B6'
b6_cor.to_csv('%s/files/protein_pairs_B6.txt' % wd, sep='\t', index=False)
print '[INFO] B6 mouse correlations done'

s9_cor = [(p1, p2, pearson(swath_quant.ix[p1, s9], swath_quant.ix[p2, s9])) for p1, p2 in subnetwork]
s9_cor = [(p1, p2, cor, pvalue, meas) for p1, p2, (cor, pvalue, meas) in s9_cor]
s9_cor = DataFrame(s9_cor, columns=['p1', 'p2', 'cor', 'pvalue', 'meas']).dropna()
s9_cor = s9_cor[s9_cor['meas'] > 2]
s9_cor['adj.pvalue'] = multipletests(s9_cor['pvalue'], method='fdr_bh')[1]
s9_cor = s9_cor.sort('adj.pvalue')
s9_cor['strain'] = 'S9'
s9_cor.to_csv('%s/files/protein_pairs_S9.txt' % wd, sep='\t', index=False)
print '[INFO] S9 mouse correlations done'

# ---- Plot pairs
p1, p2 = 'MT1_MOUSE', 'MT2_MOUSE'

plot_df = b6_cor[np.bitwise_and(b6_cor['p1'] == p1, b6_cor['p2'] == p2)]
plot_df = plot_df.append(s9_cor[np.bitwise_and(s9_cor['p1'] == p1, s9_cor['p2'] == p2)])

plot_df = swath_quant.ix[[p1, p2]]
plot_df['protein'] = plot_df.index
plot_df = melt(plot_df, id_vars='protein')
plot_df['time'] = [int(i.split('_')[1][1]) for i in plot_df['variable']]
plot_df['strain'] = [i.split('_')[0] for i in plot_df['variable']]
plot_df['condition'] = [i.split('_')[2] for i in plot_df['variable']]

grid = sns.FacetGrid(plot_df, hue='protein', col='condition', row='strain')
grid.map(plt.axhline, y=0, ls=':', c='.5')
grid.map(plt.plot, 'time', 'value', marker='', ms=4)
grid.add_legend()
plt.savefig('%s/reports/Figure8_%s_%s_pairs_correlation.pdf' % (wd, hypothesis, fdr_thres), bbox_inches='tight')
plt.close('all')