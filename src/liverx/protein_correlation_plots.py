from __future__ import division
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.lines as mlines
from sklearn.metrics.regression import mean_absolute_error
from liverx import wd
from pandas import DataFrame, Series, read_csv, melt


def rmse(x, y):
    mask = np.bitwise_and(np.isfinite(x.values), np.isfinite(y.values))
    return mean_absolute_error(x[mask], y[mask])


sns.set(style='ticks', palette='pastel', color_codes=True)

swath = read_csv('%s/data/subset-data-for-figure8.csv' % wd).append(read_csv('%s/data/result_swath_fasted_v2.3.7_fRfS_20160420.csv' % wd, index_col=0)).dropna()


# -- Defined conditions
b6_fed = ['B6_T1_fed-B6_T0_fed', 'B6_T2_fed-B6_T0_fed']
s9_fed = ['S9_T1_fed-S9_T0_fed', 'S9_T2_fed-S9_T0_fed']


# -- Plotting H2
p_pairs = [('ITB1_MOUSE', 'ICAM1_MOUSE'), ('AKC1H_MOUSE', '3BHS5_MOUSE'), ('CP4AA_MOUSE', 'UD11_MOUSE'), ('AOFB_MOUSE', 'ADH1_MOUSE')]

(f, grid), r_pos = plt.subplots(len(p_pairs), 2, figsize=(5, 3 * len(p_pairs)), sharey='row', sharex='col'), 0
for p1, p2 in p_pairs:
    plot_df = swath[[i in (p1, p2) for i in swath['Protein']]]

    plot_df[b6_fed] = plot_df[b6_fed].subtract(plot_df[b6_fed[0]], axis=0)
    plot_df[s9_fed] = plot_df[s9_fed].subtract(plot_df[s9_fed[0]], axis=0)
    plot_df[b6_fasted] = plot_df[b6_fasted].subtract(plot_df[b6_fasted[0]], axis=0)
    plot_df[s9_fasted] = plot_df[s9_fasted].subtract(plot_df[s9_fasted[0]], axis=0)

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
            ax.scatter(x, y, s=80, c=protein_colour, edgecolors='none')
            ax.axhline(y=0, lw=.3, c='gray', ls=':', alpha=.6)
            ax.set_ylabel('Fold-change')

            x, y = plot_df[plot_df.apply(lambda df: df['protein'] == protein and df['strain'] == 'S9' and df['condition'] == condition, axis=1)][['time', 'value']].T.values
            ax.plot(x, y, ls='--', c=protein_colour, label='S9')
            ax.scatter(x, y, s=80, label=protein, c=protein_colour, edgecolors='none')
            ax.axhline(y=0, lw=.3, c='gray', ls=':', alpha=.6)

            ax.set_title(condition.lower())

    sns.despine()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    r_pos += 1

plt.savefig('%s/reports/%s_%s_pairs_correlation.pdf' % (wd, hypothesis, fdr_thres), bbox_inches='tight')
plt.close('all')
print '[INFO] Gain/loss of correlation plot generated!'


# ---- Plotting H4
hypothesis, fdr_thres = ('H4', '0.05')
corr = read_csv('%s/files/protein_pairs_%s_%s.txt' % (wd, hypothesis, fdr_thres), sep='\t')
corr = corr[corr['meas'] == 6]

corr = corr[corr['e_pvalue'] < 0.05]
corr['rmse'] = [rmse(swath.ix[p1, b6_fed], swath.ix[p1, b6_fasted]) + rmse(swath.ix[p2, b6_fed], swath.ix[p2, b6_fasted]) + rmse(swath.ix[p1, s9_fed], swath.ix[p1, s9_fasted]) + rmse(swath.ix[p2, s9_fed], swath.ix[p2, s9_fasted]) for p1, p2 in corr[['p1', 'p2']].values]
corr = corr.sort('rmse', ascending=False)

(f, grid), r_pos = plt.subplots(len(corr), 3, figsize=(5, 3.5 * len(corr)), sharey='row', sharex='col'), 0
for p1, p2 in corr[['p1', 'p2']].values:

    palette = dict(zip([p1, p2], [colors.rgb2hex((r, g, b)) for r, g, b in sns.color_palette('Paired')[4:6]]))
    style = {'B6': '-', 'S9': '--'}

    for t in [0, 1, 2]:
        ax = grid[r_pos][t]
        ax.set_xticks(range(2))

        for s in ['B6', 'S9']:

            plot_df = swath.ix[[p1, p2]].eval('%s_T%d_FED - %s_T%d_FASTED' % (s, t, s, t)).to_frame(name='Fasted')
            plot_df['Fed'] = 0

            for p in [p1, p2]:

                x, y = [0, 1], list(plot_df.ix[p][['Fed', 'Fasted']])
                ax.plot(plot_df.ix[p][['Fed', 'Fasted']], ls=style[s], color=palette[p], label=s)
                ax.scatter(x, y, color=palette[p], s=80)
                ax.axhline(y=0, lw=.3, c='gray', ls=':', alpha=.6)
                ax.set_title('T%d' % t)
                ax.set_xticklabels(['Fed', 'Fasted'])
                ax.set_ylabel('Fold-change')

    sns.despine()

    handles = [mlines.Line2D([], [], color='black', linestyle=style[s], markersize=15, label=s) for s in ['B6', 'S9']]
    handles.extend(mlines.Line2D([], [], color=palette[p], linestyle='-', markersize=15, label=p) for p in [p1, p2])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), handles=handles)

    r_pos += 1

plt.savefig('%s/reports/Figure8_%s_%s_pairs_correlation.pdf' % (wd, hypothesis, fdr_thres), bbox_inches='tight')
plt.close('all')
print '[INFO] Gain/loss of correlation plot generated!'
