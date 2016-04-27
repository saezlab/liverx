import itertools as it
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from liverx import wd
from matplotlib.colors import rgb2hex
from matplotlib_venn import venn3, venn3_circles
from scipy.cluster.hierarchy import dendrogram, linkage
from pandas import DataFrame, Series, read_csv, pivot_table, melt
from sklearn.decomposition.pca import PCA
from liverx.utils import pearson

sns.set(style='ticks', palette='pastel', color_codes=True)

# Import new preprocessed data-sets
srm = read_csv('%s/data/result_srm_v2.3.7_fRfS.csv' % wd)
srm['type'] = 'SRM'
swath = read_csv('%s/data/result_swath_v2.3.7_fRfS.csv' % wd)
swath['type'] = 'SWATH'
shotgun = read_csv('%s/data/result_sg_v2.3.7_fRfS.csv' % wd)
shotgun['type'] = 'Shotgun'

srm_fc = pivot_table(srm, 'log2FC', 'Protein', 'Label', aggfunc='mean')
swath_fc = pivot_table(swath, 'log2FC', 'Protein', 'Label', aggfunc='mean')
shotgun_fc = pivot_table(shotgun, 'log2FC', 'Protein', 'Label', aggfunc='mean')

srm_quant = read_csv('%s/data/result_srm_v2.3.7_protein_quant.tab' % wd, sep='\t')
swath_quant = read_csv('%s/data/result_swath_v2.3.7_protein_quant.tab' % wd, sep='\t')
shotgun_quant = read_csv('%s/data/result_sg_v2.3.7_protein_quant.tab' % wd, sep='\t')

datasets = {'SWATH': swath_fc, 'Shotgun': shotgun_fc, 'SRM': srm_fc}
datasets_fc = {'SWATH': swath, 'Shotgun': shotgun, 'SRM': srm}
datasets_quant = {'SWATH': swath_quant, 'Shotgun': shotgun_quant, 'SRM': srm_quant}

datasets_colour = dict(zip(*(['SWATH', 'Shotgun', 'SRM'], sns.color_palette('Set1', 3))))
datasets_colour = {c: rgb2hex(datasets_colour[c]) for c in datasets_colour}

# Conditions
cond_h1 = ['B6_T1_fed-B6_T0_fed', 'B6_T2_fed-B6_T0_fed', 'B6_T2_fed-B6_T1_fed', 'S9_T1_fed-S9_T0_fed', 'S9_T2_fed-S9_T0_fed', 'S9_T2_fed-S9_T1_fed']

# -- Figure 2 - Data distributions
x_order = ['Shotgun', 'SWATH', 'SRM']

plot_df = DataFrame([(v, name) for name in x_order for v in melt(datasets[name][cond_h1])['value']], columns=['value', 'type']).dropna()
sns.violinplot(x='type', y='value', data=plot_df, order=x_order, hue_order=x_order, palette=datasets_colour)
sns.despine(trim=True)
plt.xlabel('')
plt.ylabel('Fold-change (log2)')
plt.axhline(0, c='gray', alpha=0.5, lw=.3)
plt.savefig('%s/reports/Figure2_boxplot.pdf' % wd, bbox_inches='tight')
plt.close('all')

# Distribution histograms
plot_df = [melt(datasets[name][cond_h1])['value'] for name in x_order]
plt.hist(plot_df, histtype='stepfilled', alpha=0.4, bins=80, normed=True, label=x_order, color=[datasets_colour[x] for x in x_order])
sns.despine(trim=True)
plt.legend()
plt.ylabel('Fold-change (normed log2)')
plt.savefig(wd + '/reports/Figure2_histogram.pdf', bbox_inches='tight')
plt.close('all')

# Volcanos
plot_df = [(f, -np.log10(p), l, t) for name in x_order for f, p, l, t in datasets_fc[name][['log2FC', 'adj.pvalue', 'Label', 'type']].values if l in cond_h1]
plot_df = DataFrame(plot_df, columns=['fc', 'pvalue', 'condition', 'type'])
g = sns.lmplot(x='fc', y='pvalue', data=plot_df, col='type', hue='type', col_order=x_order, palette=datasets_colour, sharex=False, sharey=False, scatter_kws={'alpha': .4})
g.set_titles(col_template='{col_name}')
g.set_axis_labels('Fold-change (log2)', 'p-value (-log10 FDR)')
sns.despine(trim=True)
plt.savefig(wd + '/reports/Figure2_volcanos.pdf', bbox_inches='tight')
plt.close('all')

# -- Overlap among data-sets
srm_set, swath_set, shotgun_set = set(srm_fc.index), set(swath_fc.index), set(shotgun_fc.index)
venn3([srm_set, swath_set, shotgun_set], set_labels=('SRM', 'SWATH', 'Shotgun'), set_colors=(datasets_colour['SRM'], datasets_colour['SWATH'], datasets_colour['Shotgun']))
plt.savefig('%s/reports/Supp_material_overlap.pdf' % wd, bbox_inches='tight')
plt.close('all')

ov_set = srm_set.intersection(swath_set).intersection(shotgun_set)

plot_df = DataFrame([(v, name) for name in x_order for v in melt(datasets[name].ix[ov_set, cond_h1])['value']], columns=['value', 'type']).dropna()
sns.violinplot(x='type', y='value', data=plot_df, order=x_order, hue_order=x_order, palette=datasets_colour)
sns.despine(trim=True)
plt.xlabel('')
plt.ylabel('Fold-change (log2)')
plt.axhline(0, c='gray', alpha=0.5, lw=.3)
plt.savefig('%s/reports/Supp_material_overlap_boxplot.pdf' % wd, bbox_inches='tight')
plt.close('all')

# Distribution histograms
plot_df = [melt(datasets[name].ix[ov_set, cond_h1])['value'] for name in x_order]
plt.hist(plot_df, histtype='stepfilled', alpha=0.4, bins=80, normed=True, label=x_order, color=[datasets_colour[x] for x in x_order])
sns.despine(trim=True)
plt.legend()
plt.ylabel('Fold-change (normed log2)')
plt.savefig(wd + '/reports/Supp_material_overlap_histogram.pdf', bbox_inches='tight')
plt.close('all')

# Volcanos
plot_df = [(f, -np.log10(p), l, t) for name in x_order for prot, f, p, l, t in datasets_fc[name][['Protein', 'log2FC', 'adj.pvalue', 'Label', 'type']].values if l in cond_h1 and prot in ov_set]
plot_df = DataFrame(plot_df, columns=['fc', 'pvalue', 'condition', 'type'])
g = sns.lmplot(x='fc', y='pvalue', data=plot_df, col='type', hue='type', col_order=x_order, palette=datasets_colour, sharex=False, sharey=False, scatter_kws={'alpha': .4})
g.set(xlim=(-6, 6))
g.set_axis_labels('Fold-change (log2)', 'p-value (-log10 FDR)')
sns.despine(trim=True)
plt.savefig(wd + '/reports/Supp_material_overlap_volcanos.pdf', bbox_inches='tight')
plt.close('all')


# -- Figure 3
cond_h1_aux = list(set(cond_h1).difference(['B6_T2_fed-B6_T1_fed', 'S9_T2_fed-S9_T1_fed']))
for name, dataset in datasets_fc.items():
    dataset_fc = pivot_table(dataset, values='log2FC', index='Protein', columns='Label')[cond_h1_aux]
    dataset_pvalue = pivot_table(dataset, values='adj.pvalue', index='Protein', columns='Label')[cond_h1_aux]

    dataset_fc = dataset_fc.ix[dataset_pvalue[(dataset_pvalue[cond_h1_aux] < 0.05).sum(1) > 0].index, cond_h1_aux]

    sns.clustermap(dataset_fc.corr())
    plt.savefig('%s/reports/Figure3_correlations_%s.pdf' % (wd, name), bbox_inches='tight')
    plt.close('all')

    print '[INFO] ' + name + ': ' + str(dataset_fc.shape)


fed = ['B6_T0_FED', 'B6_T1_FED', 'B6_T2_FED', 'S9_T0_FED', 'S9_T1_FED', 'S9_T2_FED']
for name, dataset in datasets_fc.items():
    dataset_fc = pivot_table(dataset, values='log2FC', index='Protein', columns='Label')[cond_h1_aux]
    dataset_pvalue = pivot_table(dataset, values='adj.pvalue', index='Protein', columns='Label')[cond_h1_aux]
    dataset_quant = datasets_quant[name]

    dataset_quant = dataset_quant.ix[dataset_pvalue[(dataset_pvalue[cond_h1_aux] < 0.05).sum(1) > 0].index, fed]

    sns.clustermap(dataset_quant.corr(), cmap=sns.light_palette((210, 90, 60), input='husl', as_cmap=True))
    plt.savefig('%s/reports/Figure3_correlations_protein_estimates_%s.pdf' % (wd, name), bbox_inches='tight')
    plt.close('all')

    print '[INFO] ' + name + ': ' + str(dataset_quant.shape)


# -- Figure 4
cond_h1_aux = list(set(cond_h1).difference(['B6_T2_fed-B6_T1_fed', 'S9_T2_fed-S9_T1_fed']))

sns.set(style='ticks', palette='deep')
plt.figure(figsize=(20, 15))
m_plot, pos = plt.GridSpec(3, 4, wspace=.6, hspace=.6), 0

comp_corr_df = {}
for comp1, comp2 in it.combinations(datasets_fc.keys(), 2):
    dataset_fc_comp1 = pivot_table(datasets_fc[comp1], values='log2FC', index='Protein', columns='Label')[cond_h1]
    dataset_fc_comp2 = pivot_table(datasets_fc[comp2], values='log2FC', index='Protein', columns='Label')[cond_h1]

    dataset_pvalue_comp1 = pivot_table(datasets_fc[comp1], values='adj.pvalue', index='Protein', columns='Label')[cond_h1]
    dataset_pvalue_comp2 = pivot_table(datasets_fc[comp2], values='adj.pvalue', index='Protein', columns='Label')[cond_h1]

    dataset_fc_comp1 = dataset_fc_comp1.ix[dataset_pvalue_comp1[(dataset_pvalue_comp1[cond_h1_aux] < 0.05).sum(1) > 0].index, cond_h1_aux]
    dataset_fc_comp2 = dataset_fc_comp2.ix[dataset_pvalue_comp2[(dataset_pvalue_comp2[cond_h1_aux] < 0.05).sum(1) > 0].index, cond_h1_aux]

    overlap_proteins = set(dataset_fc_comp1.index).intersection(dataset_fc_comp2.index)
    dataset_fc_comp1 = dataset_fc_comp1.ix[overlap_proteins]
    dataset_fc_comp2 = dataset_fc_comp2.ix[overlap_proteins]

    comp_corr = dataset_fc_comp1.corrwith(dataset_fc_comp2)

    comp_corr_df['%s, %s' % (comp1, comp2)] = comp_corr.to_dict()

    # Plot scatter
    for cond in cond_h1_aux:
        ax = plt.subplot(m_plot[pos])

        cor, pval, n_meas = pearson(dataset_fc_comp1.ix[overlap_proteins, cond], dataset_fc_comp2.ix[overlap_proteins, cond])

        sns.regplot(dataset_fc_comp1.ix[overlap_proteins, cond], dataset_fc_comp2.ix[overlap_proteins, cond], color='#95a5a6', ax=ax)
        sns.despine(ax=ax)
        ax.set_title('%s\nPearson: %.2f; p-val: %.2e; proteins: %d' % (cond, cor, pval, n_meas))
        ax.set_xlabel(comp1)
        ax.set_ylabel(comp2)
        ax.axhline(0, ls='--', lw=0.3, c='gray', alpha=.6)
        ax.axvline(0, ls='--', lw=0.3, c='gray', alpha=.6)

        pos += 1

plt.savefig('%s/reports/Supp_material_correlations_scatter.pdf' % wd, bbox_inches='tight')
plt.close('all')

# Plot correlation coefficients
comp_corr_df = DataFrame(comp_corr_df)
comp_corr_df['index'] = comp_corr_df.index
comp_corr_df = melt(comp_corr_df, id_vars='index')

g = sns.FacetGrid(comp_corr_df, col='variable')
g.map(plt.axvline, x=0, ls=':', c='.5', color='#95a5a6')
g.map(sns.pointplot, 'value', 'index', marker='o', ms=4, color='#95a5a6', linestyles='')
g.set_titles(col_template='{col_name}')
g.set_axis_labels('Correlation coefficient', '')
plt.savefig('%s/reports/Figure4_correlations_%s.pdf' % (wd, 'SRM'), bbox_inches='tight')
plt.close('all')


# -- Figure 5
(f, m_plot), pos = plt.subplots(3, 3, figsize=(10, 15)), 0
for name, dataset in datasets_quant.items():
    plot_df = dataset.loc[:, ['FED' in i.upper() for i in dataset.columns]].copy()

    ax = m_plot[pos][0]
    dendrogram(linkage(plot_df.T, method='complete', metric='euclidean'), ax=ax, labels=plot_df.columns, orientation='left')
    sns.despine(trim=True, ax=ax)
    ax.set_title('%s (%s)' % (name, 'euclidean'))

    ax = m_plot[pos][1]
    sns.despine(trim=True, ax=ax)
    ax.set_xticklabels('')
    ax.set_yticklabels('')

    ax = m_plot[pos][2]
    dendrogram(linkage(plot_df.T, method='complete', metric='correlation'), ax=ax, labels=plot_df.columns, orientation='left')
    sns.despine(trim=True, ax=ax)
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


# ---- Figure 6
(f, m_plot), pos = plt.subplots(3, 2, sharex=False, sharey=False, figsize=(12, 22)), 0
for name, dataset in datasets_quant.items():
    plot_df = dataset.loc[:, ['FED' in i.upper() for i in dataset.columns]].T

    n_components = 3
    pca_o = PCA(n_components=n_components).fit(plot_df)
    pcs = pca_o.transform(plot_df)
    explained_var = ['%.2f' % (pca_o.explained_variance_ratio_[i] * 100) for i in range(n_components)]

    # Plot 1
    ax = m_plot[pos][0]
    x_pc, y_pc = 0, 1
    ax.scatter(pcs[:, x_pc], pcs[:, y_pc], s=90, c=datasets_colour[name], linewidths=0)
    ax.set_xlabel('PC 1 (%s%%)' % explained_var[x_pc])
    ax.set_ylabel('PC 2 (%s%%)' % explained_var[y_pc])
    ax.set_title(name)
    sns.despine(ax=ax)

    for i, txt in enumerate(plot_df.index):
        ax.annotate(txt, (pcs[:, x_pc][i], pcs[:, y_pc][i]), size='x-small')

    # Plot 2
    ax = m_plot[pos][1]
    x_pc, y_pc = 0, 2
    ax.scatter(pcs[:, x_pc], pcs[:, y_pc], s=90, c=datasets_colour[name], linewidths=0)
    ax.set_xlabel('PC 1 (%s%%)' % explained_var[x_pc])
    ax.set_ylabel('PC 3 (%s%%)' % explained_var[y_pc])
    ax.set_title(name)
    sns.despine(ax=ax)

    for i, txt in enumerate(plot_df.index):
        ax.annotate(txt, (pcs[:, x_pc][i], pcs[:, y_pc][i]), size='x-small')

    pos += 1

plt.savefig(wd + '/reports/Figure6_pca.pdf', bbox_inches='tight')
plt.close('all')
