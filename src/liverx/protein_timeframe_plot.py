import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from liverx import wd
from matplotlib.colors import rgb2hex
from pandas import DataFrame, Series, read_csv, melt

style = {'B6': '-', 'S9': '--'}

# Conditions
b6_fed, s9_fed = ['B6_T1_fed-B6_T0_fed', 'B6_T2_fed-B6_T0_fed'], ['S9_T1_fed-S9_T0_fed', 'S9_T2_fed-S9_T0_fed']

# -- Figure 5
# Import data-sets
swath = read_csv('%s/data/subset-data-for-figure5.csv' % wd)
plist = list(read_csv('%s/files/protein_list.txt' % wd, sep='\t', header=None)[0])

# Plotting
colour = rgb2hex(sns.color_palette('Paired')[1])

sns.set(style='ticks', palette='pastel', color_codes=True)
(f, grid), r_pos = plt.subplots(len(plist) / 2, 2, figsize=(5, 2 * (len(plist) / 2)), sharey=False, sharex=True, ), 0
f.tight_layout()
f.subplots_adjust(wspace=.3, hspace=.3)
for p in plist:
    plot_df = swath[swath['Protein'] == p]

    ax = grid[r_pos / 2][r_pos % 2]
    ax.set_xticks(range(3))

    for label, conditions in [('B6', b6_fed), ('S9', s9_fed)]:
        x, y = [0, 1, 2], [0] + list(plot_df.loc[[i in conditions for i in plot_df['Label']], 'log2FC'])
        se = [0] + list(plot_df.loc[[i in conditions for i in plot_df['Label']], 'SE'])

        ax.plot(x, y, ls=style[label], c=colour, label=label)
        for xs, ys, ses in zip(x, y, se):
            ax.scatter(xs, ys, s=60, c=colour, edgecolors='white')
            ax.errorbar(xs, ys, ses, c=colour)

        ax.axhline(y=0, lw=.3, c='gray', ls='-')
        ax.set_title(p)

    sns.despine()

    if (r_pos - 2) >= (len(plist) / 2):
        ax.set_xlabel('Weeks')
        ax.set_xticklabels(['0', '6', '12'])

    if r_pos % 2 == 1:
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    else:
        ax.set_ylabel('Fold-change')

    r_pos += 1

plt.savefig('%s/reports/Figure5.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Protein time series plot generated!'


# -- Figure 7a
p_pairs = [('ITB1_MOUSE', 'ICAM1_MOUSE'), ('AKC1H_MOUSE', '3BHS5_MOUSE'), ('CP4AA_MOUSE', 'UD11_MOUSE'), ('AOFB_MOUSE', 'ADH1_MOUSE')]

sns.set(style='ticks', palette='pastel', color_codes=True)
(f, grid), r_pos = plt.subplots(len(p_pairs) / 2, 2, figsize=(6, 2 * (len(p_pairs) / 2)), sharey=False, sharex=True), 0
f.tight_layout()
f.subplots_adjust(wspace=1., hspace=.3)
for pair in p_pairs:

    ax = grid[r_pos / 2][r_pos % 2]
    ax.set_xticks(range(3))

    for p, colour in zip(*(pair, sns.color_palette('Paired').as_hex()[:2])):
        plot_df = swath[swath['Protein'] == p]

        for label, conditions in [('B6', b6_fed), ('S9', s9_fed)]:
            x, y = [0, 1, 2], [0] + list(plot_df.loc[[i in conditions for i in plot_df['Label']], 'log2FC'])
            se = [0] + list(plot_df.loc[[i in conditions for i in plot_df['Label']], 'SE'])

            ax.plot(x, y, ls=style[label], c=colour, label='%s - %s' % (label, p.split('_')[0]))
            for xs, ys, ses in zip(x, y, se):
                ax.scatter(xs, ys, s=60, c=colour, edgecolors='white')
                ax.errorbar(xs, ys, ses, c=colour)

            ax.axhline(y=0, lw=.3, c='gray', ls='-')

    ax.set_title('%s ~ %s' % (pair[0].split('_')[0], pair[1].split('_')[0]))
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    sns.despine()

    if r_pos in [2, 3]:
        ax.set_xlabel('Weeks')
        ax.set_xticklabels(['0', '6', '12'])

    if r_pos % 2 != 1:
        ax.set_ylabel('Fold-change')

    r_pos += 1

plt.savefig('%s/reports/Figure7a.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Protein time series plot generated!'


# -- Figure 7b
p_pairs = [('PAFA2_MOUSE', 'CEPT1_MOUSE'), ('SKP1_MOUSE', 'ELOC_MOUSE')]
swath = read_csv('%s/data/subset-data-for-figure8.csv' % wd).dropna()

t0, t2 = ['B6_T0_fasted - B6_T0_fed', 'S9_T0_fasted - S9_T0_fed'], ['B6_T2_fasted - B6_T2_fed', 'S9_T2_fasted - S9_T2_fed']

sns.set(style='ticks', palette='pastel', color_codes=True)
(f, grid), r_pos = plt.subplots(len(p_pairs), 2, figsize=(4, 4), sharey=True, sharex=True), 0
f.tight_layout()
f.subplots_adjust(wspace=.1, hspace=.2)
for pair in p_pairs:
    c_pos = 0

    for title, condition in zip(*(['0 weeks', '12 weeks'], [t0, t2])):
        ax = grid[r_pos][c_pos]
        ax.set_xticks(range(3))

        if r_pos == 0:
            ax.set_title(title)

        if c_pos == 0:
            ax.set_ylabel('Fold-change')

        for p, colour in zip(*(pair, sns.color_palette('Paired').as_hex()[4:6])):
            plot_df = swath[swath['Protein'] == p].set_index('Label').ix[condition]

            for strain in ['B6', 'S9']:
                x, y = [0, 1], [0] + list(plot_df.loc[[strain in i for i in plot_df.index], 'log2FC'])
                se = [0] + list(plot_df.loc[[strain in i for i in plot_df.index], 'SE'])

                ax.plot(x, y, ls=style[strain], c=colour, label='%s - %s' % (strain, p.split('_')[0]))
                for xs, ys, ses in zip(x, y, se):
                    ax.scatter(xs, ys, s=60, c=colour, edgecolors='white')
                    ax.errorbar(xs, ys, ses, c=colour)

                ax.axhline(y=0, lw=.3, c='gray', ls='-')

        if (r_pos + 1) == len(p_pairs):
            ax.set_xticklabels(['Fed', 'Fasted'])

        c_pos += 1

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    sns.despine()

    r_pos += 1

plt.savefig('%s/reports/Figure7b.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Protein time series plot generated!'


# -- Figure pathway
# Import data-sets
swath_fed = read_csv('%s/data/subset-data-for-figure5.csv' % wd).dropna()
swath_fasted = read_csv('%s/data/result_swath_fasted_v2.3.7_fRfS_20160420.csv' % wd, index_col=0).dropna()
swath = swath_fed.append(swath_fasted)

# Plotting
colour = rgb2hex(sns.color_palette('Paired')[1])

conditions = {
    'B6': {
        'fed': ['B6_T1_fed-B6_T0_fed', 'B6_T2_fed-B6_T0_fed'],
        'fasted': ['B6_T1_fasted-B6_T0_fasted', 'B6_T2_fasted-B6_T0_fasted']
    },
    'S9': {
        'fed': ['S9_T1_fed-S9_T0_fed', 'S9_T2_fed-S9_T0_fed'],
        'fasted': ['S9_T1_fasted-S9_T0_fasted', 'S9_T2_fasted-S9_T0_fasted']
    }
}


plist = list(set(read_csv('%s/files/protein_list_2.txt' % wd, sep='\t', header=None)[0]).intersection(swath['Protein']))

sns.set(style='ticks', palette='pastel', color_codes=True)
(f, grid), r_pos = plt.subplots(len(plist), 2, figsize=(4, 20), sharey=False, sharex=False), 0
f.tight_layout()
f.subplots_adjust(wspace=.6, hspace=1.)
for p in plist:
    c_pos = 0

    plot_df = swath[swath['Protein'] == p].set_index('Label')

    for condition in ['fed', 'fasted']:
        ax = grid[r_pos][c_pos]
        ax.set_xticks(range(3))

        for strain in ['B6', 'S9']:
            x, y = [0, 1, 2], [0] + list(plot_df.loc[conditions[strain][condition], 'log2FC'])
            se = [0] + list(plot_df.loc[conditions[strain][condition], 'SE'])

            ax.plot(x, y, ls=style[strain], c=colour, label=strain)
            for xs, ys, ses in zip(x, y, se):
                ax.scatter(xs, ys, s=60, c=colour, edgecolors='white')
                ax.errorbar(xs, ys, ses, c=colour)

            ax.axhline(y=0, lw=.3, c='gray', ls='-')

        ax.set_title(condition.capitalize())
        sns.despine()

        ax.set_xlabel('Weeks')
        ax.set_xticklabels(['0', '6', '12'])

        if c_pos == 1:
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        else:
            ax.set_ylabel('%s\nfold-change' % p.split('_')[0])

        c_pos += 1

    r_pos += 1

plt.gcf().set_size_inches(4, 80)
plt.savefig('%s/reports/Figure_pathway.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Protein time series plot generated!'
