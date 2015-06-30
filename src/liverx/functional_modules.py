from __future__ import division
import igraph
import itertools as it
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from sklearn.metrics.pairwise import pairwise_distances
from statsmodels.stats.multitest import multipletests
from pandas.stats.misc import zscore
from sklearn.cluster import AgglomerativeClustering, AffinityPropagation
from liverx import wd
from scipy.stats.stats import pearsonr
from pandas import DataFrame, read_csv, melt


def pearson(x, y):
    mask = np.bitwise_and(np.isfinite(x), np.isfinite(y))
    cor, pvalue = pearsonr(x[mask], y[mask])
    return cor, pvalue, sum(mask)


sns.set_style('white')

swath_quant = read_csv('%s/data/result_swath_v2.3.7_protein_quant.tab' % wd, sep='\t').replace(0.0, np.NaN)

# ---- Defined conditions
b6 = [c for c in swath_quant.columns if c.startswith('B6')]
s9 = [c for c in swath_quant.columns if c.startswith('S9')]

# ---- Import String network
network = read_csv('%s/files/string_mouse_network_filtered_800.txt' % wd, sep='\t', names=['p1', 'p2'], skiprows=1)

network_i = igraph.Graph(directed=False)
network_i.add_vertices(list(set(network['p1']).union(network['p2'])))
network_i.add_edges([(p1, p2) for p1, p2 in set(zip(*(network['p1'], network['p2'])))])
print '[INFO] String network: ', network_i.summary()

network_i = network_i.subgraph(set(network_i.vs['name']).intersection(swath_quant.index))
print '[INFO] Swath measured simplified string network: ', network_i.summary()

hypothesis, fdr_thres = 'H2', '0.05'

# ---- Import subnetwork
subnetwork = read_csv('%s/files/network_enrichment/%s_%s_network.sif' % (wd, hypothesis, fdr_thres), header=None, sep='\t', names=['p1', 'i', 'p2'])
subnetwork_i = network_i.subgraph(set(it.chain(*network_i.neighborhood(set(subnetwork['p1']).union(subnetwork['p2']), order=0))))
print '[INFO] Swath measured simplified string subnetwork: ', subnetwork_i.summary()

plot_df = swath_quant.ix[subnetwork_i.vs['name'], b6].dropna()

aff_prop = AffinityPropagation()
aff_prop.fit(plot_df)
len(aff_prop.cluster_centers_)

model = AgglomerativeClustering(n_clusters=9, linkage='average', affinity='l1')
label = model.fit_predict(plot_df)
print label

(f, grid), pos = plt.subplots(9, 1, figsize=(6, 4 * 9)), 0
for i in range(9):
    plot_df.ix[plot_df.index[label == i]].T.plot(legend=False, ax=grid[pos])
    pos += 1
plt.savefig('%s/reports/functional_modules_%s_%s.pdf' % (wd, hypothesis, fdr_thres), bbox_inches='tight')
plt.close('all')

plot_df = DataFrame([(p, pairwise_distances(swath_quant.ix[p, b6], swath_quant.ix[p, s9], metric='cosine')[0][0]) for p in plot_df.index])

plot_df = swath_quant.ix[subnetwork_i.vs['name']].dropna()
sns.clustermap(plot_df, xticklabels=True, col_cluster=False, yticklabels=False, cmap='Blues')

plt.savefig('%s/reports/functional_modules_%s_%s.pdf' % (wd, hypothesis, fdr_thres), bbox_inches='tight')
plt.close('all')
print '[INFO] Plotted'