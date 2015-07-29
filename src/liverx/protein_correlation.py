from __future__ import division
import igraph
import itertools as it
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.stats.stats import pearsonr
from liverx import wd
from liverx.utils import pearson
from statsmodels.stats.multitest import multipletests
from pandas import DataFrame, read_csv, melt


sns.set_style('ticks')

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


# ---- Randomly suffle kinase matrix while keeping NaNs position
def randomise_matrix(matrix):
    random_df = matrix.copy()
    movers = ~np.isnan(random_df.values)
    random_df.values[movers] = np.random.permutation(random_df.values[movers])
    return random_df

n_permutations = 1000
swath_rand_b6 = [randomise_matrix(swath_quant[b6]) for i in xrange(n_permutations)]
swath_rand_s9 = [randomise_matrix(swath_quant[s9]) for i in xrange(n_permutations)]
print '[INFO] SWATH data-set randomisation done: ', len(swath_rand_s9)

for hypothesis, fdr_thres in [('H2', '0.05'), ('H4', '0.05')]:
    print '[INFO] Hypothesis, FDR: ', hypothesis, fdr_thres

    # ---- Import subnetwork
    subnetwork = read_csv('%s/files/network_enrichment/%s_%s_network.sif' % (wd, hypothesis, fdr_thres), header=None, sep='\t', names=['p1', 'i', 'p2'])
    subnetwork_i = network_i.subgraph(set(it.chain(*network_i.neighborhood(set(subnetwork['p1']).union(subnetwork['p2']), order=0)))).simplify()
    print '[INFO] Swath measured simplified string subnetwork: ', subnetwork_i.summary()

    subnetwork = [tuple(subnetwork_i.vs[[e.source, e.target]]['name']) for e in subnetwork_i.es]

    def correlation(p1, p2):
        print '[INFO] %s - %s' % (p1, p2)

        # Calculate protein 1 and protein 2 correlation in both strains
        cor_b6, pvalue_b5, meas_b6 = pearson(swath_quant.ix[p1, b6], swath_quant.ix[p2, b6])
        cor_s9, pvalue_s9, meas_s9 = pearson(swath_quant.ix[p1, s9], swath_quant.ix[p2, s9])

        n_meas = min(meas_b6, meas_s9)

        if n_meas > 4:
            # Calculate correlation difference
            cor_diff = abs(cor_b6 - cor_s9)

            # Calculate random correlation difference
            rand_cor_diff = [abs(pearson(swath_rand_b6[i].ix[p1, b6], swath_rand_b6[i].ix[p2, b6])[0] - pearson(swath_rand_s9[i].ix[p1, s9], swath_rand_s9[i].ix[p2, s9])[0]) for i in xrange(n_permutations)]

            # Calculate correlation difference empirical p-value
            count = sum([(r_cor >= cor_diff >= 0) or (r_cor <= cor_diff < 0) for r_cor in rand_cor_diff])
            e_pvalue = 1 / n_permutations if count == 0 else count / n_permutations

            print '\t', p1, p2, cor_diff, e_pvalue, n_meas
            return p1, p2, cor_diff, e_pvalue, n_meas

        else:
            print '\t', p1, p2, 'NaN', 'NaN', n_meas
            return p1, p2, np.NaN, np.NaN, n_meas

    # ---- Compute correlations differences
    correlation_df = [correlation(p1, p2) for p1, p2 in subnetwork]
    correlation_df = DataFrame(correlation_df, columns=['p1', 'p2', 'diff', 'e_pvalue', 'meas']).dropna()
    correlation_df['adj_e_pvalue'] = multipletests(correlation_df['e_pvalue'], method='fdr_bh')[1]
    correlation_df = correlation_df.sort('diff', ascending=False)
    correlation_df.to_csv('%s/files/protein_pairs_%s_%s.txt' % (wd, hypothesis, fdr_thres), index=False, sep='\t')
    print '[INFO] Correlation differences calculated: ', len(correlation_df)
