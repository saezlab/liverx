import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from liverx import wd
from matplotlib.colors import rgb2hex
from pandas import DataFrame, Series, read_csv, melt


# ---- Import data-sets
swath = read_csv('%s/data/result_swath_v2.3.7_protein_quant.tab' % wd, sep='\t').replace(0.0, np.NaN)
plist = list(set(read_csv('%s/files/protein_list.txt' % wd, sep='\t', header=None)[0]).intersection(swath.index))

# ---- Defined conditions
b6 = [c for c in swath.columns if c.startswith('B6')]
s9 = [c for c in swath.columns if c.startswith('S9')]

b6_fed, b6_fasted = [c for c in b6 if c.endswith('FED')], [c for c in b6 if c.endswith('FASTED')]
s9_fed, s9_fasted = [c for c in s9 if c.endswith('FED')], [c for c in s9 if c.endswith('FASTED')]

# ---- Plotting
colour = rgb2hex(sns.color_palette('Paired')[1])
style = {'B6': '-', 'S9': '--'}

sns.set(style='ticks', palette='pastel', color_codes=True)
(f, grid), r_pos = plt.subplots(len(plist), 1, figsize=(3, 3 * len(plist)), sharey=False, sharex=True), 0
for p in plist:
    plot_df = swath.ix[p]

    plot_df[b6_fed] = plot_df[b6_fed].subtract(plot_df[b6_fed[0]])
    plot_df[s9_fed] = plot_df[s9_fed].subtract(plot_df[s9_fed[0]])

    ax = grid[r_pos]
    ax.set_xticks(range(3))

    for label, conditions in [('B6', b6_fed), ('S9', s9_fed)]:

        x, y = [0, 1, 2], plot_df[conditions]
        ax.plot(x, y, ls=style[label], c=colour, label=label)
        ax.scatter(x, y, s=80, c=colour, edgecolors='none')
        ax.axhline(y=0, lw=.3, c='gray', ls=':', alpha=.6)
        ax.set_ylabel('Fold-change')
        ax.set_title(p)

    sns.despine()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    r_pos += 1

plt.savefig('%s/reports/Protein_time_series.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Protein time series plot generated!'
