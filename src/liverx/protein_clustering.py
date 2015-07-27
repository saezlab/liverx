import igraph
import numpy as np
import itertools as it
import seaborn as sns
import matplotlib.pyplot as plt
from liverx import wd
from pandas import DataFrame, Series, read_csv

swath_quant = read_csv('%s/data/result_swath_v2.3.7_protein_quant.tab' % wd, sep='\t').replace(0.0, np.NaN)

# ---- Defined conditions
b6 = [c for c in swath_quant.columns if c.startswith('B6')]
s9 = [c for c in swath_quant.columns if c.startswith('S9')]

hypothesis, fdr_thres = ('H2', '0.05')
print '[INFO] Hypothesis, FDR: ', hypothesis, fdr_thres

# ---- Import String network
network = read_csv('%s/files/string_mouse_network_filtered_800.txt' % wd, sep='\t', names=['p1', 'p2'], skiprows=1)

network_i = igraph.Graph(directed=False)
network_i.add_vertices(list(set(network['p1']).union(network['p2'])))
network_i.add_edges([(p1, p2) for p1, p2 in set(zip(*(network['p1'], network['p2'])))])
print '[INFO] String network: ', network_i.summary()

network_i = network_i.subgraph(set(network_i.vs['name']).intersection(swath_quant.index))
print '[INFO] Swath measured simplified string network: ', network_i.summary()

# ---- Import subnetwork
subnetwork = read_csv('%s/files/network_enrichment/%s_%s_network.sif' % (wd, hypothesis, fdr_thres), header=None, sep='\t', names=['p1', 'i', 'p2'])
subnetwork_i = network_i.subgraph(set(it.chain(*network_i.neighborhood(set(subnetwork['p1']).union(subnetwork['p2']), order=0)))).simplify()
print '[INFO] Swath measured simplified string subnetwork: ', subnetwork_i.summary()

b6_fed = [c for c in b6 if c.endswith('FED')]
plot_df = swath_quant.ix[list(subnetwork_i.vs['name']), b6_fed]
plot_df = np.log2(plot_df.div(plot_df['B6_T0_FED'], axis=0))

sns.clustermap(plot_df.dropna(), col_cluster=False, yticklabels=False)
