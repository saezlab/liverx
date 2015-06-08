import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr, spearmanr
from pandas import DataFrame, Series, read_csv, pivot_table
from statsmodels.stats.multitest import multipletests
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

data, wd = '/Users/emanuel/Projects/data/liverx_mouse', '/Users/emanuel/Projects/projects/liverx'

sns.set_style('white')

# Import mapping between old and new conditions
conditions_map = read_csv(wd + '/files/conditions_map.tab', sep='\t', index_col=0).to_dict()['old']

# Import new preprocessed data-sets
srm = read_csv('%s/result_srm_v2.3.7_fRfS.csv' % data)
swath = read_csv('%s/result_swath_v2.3.7_fRfS.csv' % data)
shotgun = read_csv('%s/result_sg_v2.3.7_fRfS.csv' % data)

srm_fc = pivot_table(srm, 'log2FC', 'Protein', 'Label')
swath_fc = pivot_table(swath, 'log2FC', 'Protein', 'Label')
shotgun_fc = pivot_table(shotgun, 'log2FC', 'Protein', 'Label')

# Import shotgun old conditions
old_shotgun_fc = DataFrame([read_csv('%s/oldshotgundataset/%s' % (data, old), index_col=1)['logFC'] for new, old in conditions_map.items()]).T
old_shotgun_fc.columns = conditions_map.keys()

old_shotgun_fc_pvalue = DataFrame([read_csv('%s/oldshotgundataset/%s' % (data, old), index_col=1)['adj.pvalue'] for new, old in conditions_map.items()]).T
old_shotgun_fc_pvalue.columns = conditions_map.keys()
old_shotgun_fc_pvalue = old_shotgun_fc_pvalue.ix[old_shotgun_fc.index, old_shotgun_fc_pvalue.columns]

# Correlation old vs new
proteins = old_shotgun_fc.index
for k in conditions_map:
    grid = sns.JointGrid(shotgun_fc.ix[proteins, k].values, old_shotgun_fc.ix[proteins, k].values, ratio=50)
    grid.plot_joint(plt.scatter, color='#95a5a6')
    grid.annotate(pearsonr, loc=4)
    grid.plot_marginals(sns.rugplot, color='#95a5a6')

    grid.set_axis_labels('new shotgun', 'old shotgun')

    plt.savefig('%s/reports/new_vs_old_shotgun_%s.pdf' % (wd, k), bbox_inches='tight')
    plt.close('all')

# Correlation with adj-pvalue filters
pfilter_df = []
for pfilter in np.arange(0, 1, 0.01):
    for k in conditions_map:
        y = old_shotgun_fc[old_shotgun_fc_pvalue < pfilter][k].dropna()
        x = shotgun_fc.ix[y.index, k]

        cor, pvalue = pearsonr(x, y)

        pfilter_df.append((k, str(pfilter), cor, pvalue))

pfilter_df = DataFrame(pfilter_df, columns=['condition', 'filter', 'correlation', 'pvalue'])
pfilter_df['adj_pvalue'] = multipletests(pfilter_df['pvalue'], method='fdr_bh')[1]

g = sns.factorplot('filter', 'correlation', 'condition', pfilter_df, kind='point')
g.despine(trim=True, bottom=True)
g.set_axis_labels('p-value filter', 'correlation (pearson)')
g.set_xticklabels(step=15)

plt.savefig('%s/reports/new_vs_old_shotgun_boxplot.pdf' % wd, bbox_inches='tight')
plt.close('all')