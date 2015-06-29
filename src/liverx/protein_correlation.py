from __future__ import division
import itertools as it
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from statsmodels.stats.multitest import multipletests
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

hypothesis, fdr_thres = 'H2', '0.01'

# # ---- Import String network
# network = read_csv('%s/files/string_mouse_network_filtered_800.txt' % wd, sep='\t', names=['p1', 'p2'], skiprows=1)
# network = set(zip(*(network['p1'], network['p2'])))

# ---- Import subnetwork
subnetwork = read_csv('%s/files/network_enrichment/%s_%s_network.sif' % (wd, hypothesis, fdr_thres), header=None, sep='\t', names=['p1', 'i', 'p2'])
subnetwork = set(zip(*(subnetwork['p1'], subnetwork['p2'])))

# subnetwork = set(subnetwork['p1']).intersection(subnetwork['p2'])
# subnetwork = {(p1, p2) for p1, p2 in network if p1 in subnetwork or p2 in subnetwork}
# subnetwork = {(p1, p2) for p1, p2 in network if p1 in swath_quant.index and p2 in swath_quant.index}


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


# def empirical_pvalue(p1, p2, conditions, randomised_swath):
#     # Calculate analytical correlaiton p-value
#     cor, a_pvalue, n_meas = pearson(swath_quant.ix[p1, conditions], swath_quant.ix[p2, conditions])
#
#     # Calculate correlaiton of randomised matrices
#     random_cor = [pearson(r_matrix.ix[p1, conditions], r_matrix.ix[p2, conditions])[0] for r_matrix in randomised_swath]
#     count = sum([(r_cor >= cor >= 0) or (r_cor <= cor < 0) for r_cor in random_cor])
#
#     # Calculate empirical p-value
#     e_pvalue = 1 / n_permutations if count == 0 else count / n_permutations
#
#     print '[INFO] ', p1, p2, cor, a_pvalue, e_pvalue, n_meas
#
#     return p1, p2, cor, a_pvalue, e_pvalue, n_meas


def correlation(p1, p2):
    # Calculate protein 1 and protein 2 correlation in both strains
    cor_b6, pvalue_b5, meas_b6 = pearson(swath_quant.ix[p1, b6], swath_quant.ix[p2, b6])
    cor_s9, pvalue_s9, meas_s9 = pearson(swath_quant.ix[p1, s9], swath_quant.ix[p2, s9])

    # Calculate correlation difference
    cor_diff = abs(cor_b6 - cor_s9)

    # Calculate random correlation difference
    rand_cor_diff = [abs(pearson(swath_rand_b6[i].ix[p1, b6], swath_rand_b6[i].ix[p2, b6])[0] - pearson(swath_rand_s9[i].ix[p1, s9], swath_rand_s9[i].ix[p2, s9])[0]) for i in xrange(n_permutations)]

    # Calculate correlation difference empirical p-value
    count = sum([(r_cor >= cor_diff >= 0) or (r_cor <= cor_diff < 0) for r_cor in rand_cor_diff])
    e_pvalue = 1 / n_permutations if count == 0 else count / n_permutations

    print '[INFO] ', p1, p2, cor_diff, e_pvalue, min(meas_b6, meas_s9)

    return p1, p2, cor_diff, e_pvalue, min(meas_b6, meas_s9)

# ---- Compute correlations differences
# b6_cor = [empirical_pvalue(p1, p2, b6, swath_rand_b6) for p1, p2 in subnetwork]
# b6_cor = DataFrame(b6_cor, columns=['p1', 'p2', 'cor', 'a_pvalue', 'e_pvalue', 'meas']).dropna()
# b6_cor = b6_cor[b6_cor['meas'] > 2]
# b6_cor['adj_a_pvalue'] = multipletests(b6_cor['a_pvalue'], method='fdr_bh')[1]
# b6_cor['adj_e_pvalue'] = multipletests(b6_cor['e_pvalue'], method='fdr_bh')[1]
# b6_cor.to_csv('%s/files/protein_pairs_B6.txt' % wd, sep='\t', index=False)
# print '[INFO] B6 mouse correlations done'
#
# s9_cor = [empirical_pvalue(p1, p2, s9, swath_rand_s9) for p1, p2 in subnetwork]
# s9_cor = DataFrame(s9_cor, columns=['p1', 'p2', 'cor', 'a_pvalue', 'e_pvalue', 'meas']).dropna()
# s9_cor = s9_cor[s9_cor['meas'] > 2]
# s9_cor['adj_a_pvalue'] = multipletests(s9_cor['a_pvalue'], method='fdr_bh')[1]
# s9_cor['adj_e_pvalue'] = multipletests(s9_cor['e_pvalue'], method='fdr_bh')[1]
# s9_cor.to_csv('%s/files/protein_pairs_S9.txt' % wd, sep='\t', index=False)
# print '[INFO] S9 mouse correlations done'

correlation_df = [correlation(p1, p2) for p1, p2 in subnetwork]
correlation_df = DataFrame(correlation_df, columns=['p1', 'p2', 'diff', 'e_pvalue', 'meas']).dropna()
correlation_df = correlation_df[correlation_df['meas'] > 3]
correlation_df['adj_e_pvalue'] = multipletests(correlation_df['e_pvalue'], method='fdr_bh')[1]
correlation_df = correlation_df.sort('diff', ascending=False)
print '[INFO] Correlation differences calculated: ', len(correlation_df)

# ---- Plot pairs
dif_pairs = set(zip(*correlation_df[correlation_df['e_pvalue'] < 0.01][['p1', 'p2']].T.values))

(f, grid), r_pos = plt.subplots(len(dif_pairs), 2, figsize=(7, 4 * len(dif_pairs)), sharey=True, sharex=True), 0
for p1, p2 in dif_pairs:

    plot_df = b6_cor[np.bitwise_and(b6_cor['p1'] == p1, b6_cor['p2'] == p2)]
    plot_df = plot_df.append(s9_cor[np.bitwise_and(s9_cor['p1'] == p1, s9_cor['p2'] == p2)])

    plot_df = swath_quant.ix[[p1, p2]]
    plot_df['protein'] = plot_df.index
    plot_df = melt(plot_df, id_vars='protein')
    plot_df['time'] = [int(i.split('_')[1][1]) for i in plot_df['variable']]
    plot_df['strain'] = [i.split('_')[0] for i in plot_df['variable']]
    plot_df['condition'] = [i.split('_')[2] for i in plot_df['variable']]

    p1_color, p2_color = [colors.rgb2hex(c) for c in sns.color_palette('Paired')[:2]]

    for c_pos, condition in [(0, 'FED'), (1, 'FASTED')]:
        ax = grid[r_pos][c_pos]
        ax.set_xticks(range(3))

        for protein, protein_colour in [(p1, p1_color), (p2, p2_color)]:
            x, y = plot_df[plot_df.apply(lambda df: df['protein'] == protein and df['strain'] == 'B6' and df['condition'] == condition, axis=1)][['time', 'value']].T.values
            ax.plot(x, y, ls='-', c=protein_colour, label='B6')
            ax.scatter(x, y, s=50, c=protein_colour, edgecolors='none')

            x, y = plot_df[plot_df.apply(lambda df: df['protein'] == protein and df['strain'] == 'S9' and df['condition'] == condition, axis=1)][['time', 'value']].T.values
            ax.plot(x, y, ls='--', c=protein_colour, label='S9')
            ax.scatter(x, y, s=50, label=protein, c=protein_colour, edgecolors='none')

            ax.set_title(condition.lower())

    sns.despine()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    r_pos += 1

plt.savefig('%s/reports/Figure8_%s_%s_pairs_correlation.pdf' % (wd, hypothesis, fdr_thres), bbox_inches='tight')
plt.close('all')