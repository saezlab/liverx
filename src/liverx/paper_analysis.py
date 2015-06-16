import itertools as it
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from liverx import wd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats.stats import pearsonr, spearmanr
from pandas import DataFrame, Series, read_csv, pivot_table, melt
from sklearn.decomposition.pca import PCA


def pearson(x, y):
    mask = np.bitwise_and(np.isfinite(x), np.isfinite(y))
    cor, pvalue = pearsonr(x[mask], y[mask]) if np.sum(mask) > 1 else (np.NaN, np.NaN)
    return cor, pvalue

sns.set_style('white')

# Import mapping between old and new conditions
condition_name = read_csv(wd + '/files/conditions_name_map.tab', sep='\t', index_col=0)

# Import new preprocessed data-sets
srm = read_csv('%s/data/result_srm_v2.3.7_fRfS.csv' % wd)
swath = read_csv('%s/data/result_swath_v2.3.7_fRfS.csv' % wd)
shotgun = read_csv('%s/data/result_sg_v2.3.7_fRfS.csv' % wd)

srm_fc = pivot_table(srm, 'log2FC', 'Protein', 'Label', aggfunc='mean')
swath_fc = pivot_table(swath, 'log2FC', 'Protein', 'Label', aggfunc='mean')
shotgun_fc = pivot_table(shotgun, 'log2FC', 'Protein', 'Label', aggfunc='mean')

srm_quant = read_csv('%s/data/result_srm_v2.3.7_protein_quant.tab' % wd, sep='\t')
swath_quant = read_csv('%s/data/result_swath_v2.3.7_protein_quant.tab' % wd, sep='\t')
shotgun_quant = read_csv('%s/data/result_sg_v2.3.7_protein_quant.tab' % wd, sep='\t')

datasets = {'SWATH': swath_fc, 'Shotgun': shotgun_fc, 'SRM': srm_fc}
datasets_fc = {'SWATH': swath, 'Shotgun': shotgun, 'SRM': srm}
datasets_quant = {'SWATH': swath_quant, 'Shotgun': shotgun_quant, 'SRM': srm_quant}

datasets_colour = {'SWATH': '#9E0B0F', 'Shotgun': '#808080', 'SRM': '#2748AB'}

# Conditions
cond_h1 = ['B6_T1_fed-B6_T0_fed', 'B6_T2_fed-B6_T0_fed', 'B6_T2_fed-B6_T1_fed', 'S9_T1_fed-S9_T0_fed', 'S9_T2_fed-S9_T0_fed', 'S9_T2_fed-S9_T1_fed']

# ---- Figure 2
# Distribution boxplots
(f, m_plot), pos = plt.subplots(1, 3, sharey=True, figsize=(5, 20)), 0
for name, dataset in datasets.items():
    ax = m_plot[pos]

    plot_df = melt(dataset[cond_h1])['value']

    sns.boxplot(plot_df, color=datasets_colour[name], ax=ax)
    sns.despine(left=True, bottom=True, ax=ax)
    ax.set_title(name)

    if pos == 0:
        ax.set_ylabel('Fold-change (log2)')

    pos += 1

plt.savefig(wd + '/reports/Figure2_boxplot.pdf', bbox_inches='tight')
plt.close('all')

# Distribution histograms
plot_df, colours, labels = zip(*[(melt(dataset[cond_h1]).dropna()['value'].values, datasets_colour[name], name) for name, dataset in datasets.items()])

plt.hist(list(plot_df), histtype='stepfilled', alpha=0.4, bins=80, normed=True, label=labels)
sns.despine(left=True, bottom=True)
plt.legend()
plt.ylabel('Fold-change (normed log2)')

plt.savefig(wd + '/reports/Figure2_histogram.pdf', bbox_inches='tight')
plt.close('all')

# Volcanos
(f, m_plot), pos = plt.subplots(1, 3, sharey=True, figsize=(15, 5)), 0
for name, dataset in datasets_fc.items():
    ax = m_plot[pos]

    plot_df = dataset[[c in cond_h1 for c in dataset['Label']]]

    x, y = plot_df['log2FC'], -np.log10(plot_df['adj.pvalue'])

    ax.scatter(x, y, linewidths=0, s=20, c=datasets_colour[name], alpha=0.5)
    sns.despine(left=True, bottom=True, ax=ax)

    ax.axhline(-np.log10(0.05), ls='--', c='#bbbbbb')
    ax.set_title(name)

    if pos == 0:
        ax.set_ylabel('Adj. p-value (-log10)')

    pos += 1

plt.savefig(wd + '/reports/Figure2_volcanos.pdf', bbox_inches='tight')
plt.close('all')

# ---- Figure 3
cond_h1_aux = list(set(cond_h1).difference(['B6_T2_fed-B6_T1_fed', 'S9_T2_fed-S9_T1_fed']))
for name, dataset in datasets_fc.items():
    dataset_fc = pivot_table(dataset, values='log2FC', index='Protein', columns='Label')[cond_h1_aux]
    dataset_pvalue = pivot_table(dataset, values='adj.pvalue', index='Protein', columns='Label')[cond_h1_aux]

    dataset_fc = dataset_fc.ix[dataset_pvalue[(dataset_pvalue < 0.05).sum(1) > 0].index]

    sns.clustermap(dataset_fc.corr())
    plt.savefig(wd + '/reports/Figure3_correlations_%s.pdf' % name, bbox_inches='tight')
    plt.close('all')

    print '[INFO] ' + name + ': ' + str(dataset_fc.shape)

# ---- Figure 4
comp_corr_df = {}
for comp1, comp2 in it.combinations(datasets_fc.keys(), 2):
    dataset_fc_comp1 = pivot_table(datasets_fc[comp1], values='log2FC', index='Protein', columns='Label')[cond_h1]
    dataset_fc_comp2 = pivot_table(datasets_fc[comp2], values='log2FC', index='Protein', columns='Label')[cond_h1]

    dataset_pvalue_comp1 = pivot_table(datasets_fc[comp1], values='adj.pvalue', index='Protein', columns='Label')[cond_h1]
    dataset_pvalue_comp2 = pivot_table(datasets_fc[comp2], values='adj.pvalue', index='Protein', columns='Label')[cond_h1]

    dataset_fc_comp1 = dataset_fc_comp1.ix[dataset_pvalue_comp1[(dataset_pvalue_comp1 < 0.05).sum(1) > 0].index]
    dataset_fc_comp2 = dataset_fc_comp2.ix[dataset_pvalue_comp2[(dataset_pvalue_comp2 < 0.05).sum(1) > 0].index]

    overlap_proteins = set(dataset_fc_comp1.index).intersection(dataset_fc_comp2.index)
    dataset_fc_comp1 = dataset_fc_comp1.ix[overlap_proteins]
    dataset_fc_comp2 = dataset_fc_comp2.ix[overlap_proteins]

    comp_corr = dataset_fc_comp1.corrwith(dataset_fc_comp2)

    comp_corr_df['%s, %s' % (comp1, comp2)] = comp_corr.to_dict()

comp_corr_df = DataFrame(comp_corr_df)
comp_corr_df['index'] = comp_corr_df.index
comp_corr_df = melt(comp_corr_df, id_vars='index')

g = sns.factorplot('index', 'value', data=comp_corr_df, row='variable', kind='point', aspect=3, linestyles='-')
sns.despine(bottom=True)
g.set_xlabels('')
g.set_ylabels('Correlation coefficient')
plt.savefig(wd + '/reports/Figure4_correlations_%s.pdf' % name, bbox_inches='tight')
plt.close('all')

# ---- Figure 5
(f, m_plot), pos = plt.subplots(3, 3, figsize=(10, 15)), 0
for name, dataset in datasets_quant.items():
    plot_df = dataset.loc[:, ['FED' in i.upper() for i in dataset.columns]].copy()

    ax = m_plot[pos][0]
    dendrogram(linkage(plot_df.T, method='complete', metric='euclidean'), ax=ax, labels=plot_df.columns, orientation='left')
    sns.despine(left=True, bottom=True, ax=ax)
    ax.set_title('%s (%s)' % (name, 'euclidean'))

    ax = m_plot[pos][1]
    sns.despine(left=True, bottom=True, ax=ax)
    ax.set_xticklabels('')
    ax.set_yticklabels('')

    ax = m_plot[pos][2]
    dendrogram(linkage(plot_df.T, method='complete', metric='correlation'), ax=ax, labels=plot_df.columns, orientation='left')
    sns.despine(left=True, bottom=True, ax=ax)
    ax.set_title('%s (%s)' % (name, 'correlation'))

    pos += 1

plt.savefig(wd + '/reports/Figure5_dendrograms.pdf', bbox_inches='tight')
plt.close('all')

cond_h1_aux = list(set(cond_h1).difference(['B6_T2_fed-B6_T1_fed', 'S9_T2_fed-S9_T1_fed']))
(f, m_plot), pos = plt.subplots(3, 3, figsize=(10, 15)), 0
for name, dataset in datasets_quant.items():
    cond = [i for i in dataset.columns if 'FED' in i.upper()]

    dataset_pvalue = pivot_table(datasets_fc[name], values='adj.pvalue', index='Protein', columns='Label')[cond_h1_aux]
    proteins = dataset_pvalue[(dataset_pvalue < 0.05).sum(1) > 0].index

    plot_df = dataset.loc[proteins, cond]

    ax = m_plot[pos][0]
    dendrogram(linkage(plot_df.T, method='complete', metric='euclidean'), ax=ax, labels=plot_df.columns, orientation='left')
    sns.despine(left=True, bottom=True, ax=ax)
    ax.set_title('%s (%s)' % (name, 'euclidean'))

    ax = m_plot[pos][1]
    sns.despine(left=True, bottom=True, ax=ax)
    ax.set_xticklabels('')
    ax.set_yticklabels('')

    ax = m_plot[pos][2]
    dendrogram(linkage(plot_df.T, method='complete', metric='correlation'), ax=ax, labels=plot_df.columns, orientation='left')
    sns.despine(left=True, bottom=True, ax=ax)
    ax.set_title('%s (%s)' % (name, 'correlation'))

    pos += 1

plt.savefig(wd + '/reports/Figure5_dendrograms_signif_protein.pdf', bbox_inches='tight')
plt.close('all')

cond_h1_aux = list(set(cond_h1).difference(['B6_T2_fed-B6_T1_fed', 'S9_T2_fed-S9_T1_fed']))
for name, dataset in datasets_quant.items():
    cond = [i for i in dataset.columns if 'FED' in i.upper()]

    dataset_pvalue = pivot_table(datasets_fc[name], values='adj.pvalue', index='Protein', columns='Label')[cond_h1_aux]
    proteins = dataset_pvalue[(dataset_pvalue < 0.05).sum(1) > 0].index

    plot_df = dataset.loc[proteins, cond]

    sns.clustermap(plot_df.corr())
    plt.savefig(wd + '/reports/Figure3_correlations_protein_estimates_%s.pdf' % name, bbox_inches='tight')
    plt.close('all')


# ---- Figure 6
(f, m_plot), pos = plt.subplots(3, 2, sharex=False, sharey=False, figsize=(12, 20)), 0
for name, dataset in datasets_quant.items():
    plot_df = dataset.loc[:, ['FED' in i.upper() for i in dataset.columns]].T

    n_components = 3
    pca_o = PCA(n_components=n_components).fit(plot_df)
    pcs = pca_o.transform(plot_df)
    explained_var = ['%.2f' % (pca_o.explained_variance_ratio_[i] * 100) for i in range(n_components)]

    # Plot 1
    ax = m_plot[pos][0]
    x_pc, y_pc = 0, 1
    ax.scatter(pcs[:, x_pc], pcs[:, y_pc], s=30, c='#95a5a6', linewidths=0)
    ax.set_xlabel('PC 1 (%s%%)' % explained_var[x_pc])
    ax.set_ylabel('PC 2 (%s%%)' % explained_var[y_pc])
    ax.set_title(name)
    sns.despine(ax=ax)

    for i, txt in enumerate(plot_df.index):
        ax.annotate(txt, (pcs[:, x_pc][i], pcs[:, y_pc][i]), size='x-small')

    # Plot 2
    ax = m_plot[pos][1]
    x_pc, y_pc = 0, 2
    ax.scatter(pcs[:, x_pc], pcs[:, y_pc], s=30, c='#95a5a6', linewidths=0)
    ax.set_xlabel('PC 1 (%s%%)' % explained_var[x_pc])
    ax.set_ylabel('PC 3 (%s%%)' % explained_var[y_pc])
    ax.set_title(name)
    sns.despine(ax=ax)

    for i, txt in enumerate(plot_df.index):
        ax.annotate(txt, (pcs[:, x_pc][i], pcs[:, y_pc][i]), size='x-small')

    pos += 1

plt.savefig(wd + '/reports/Figure6_pca.pdf', bbox_inches='tight')
plt.close('all')