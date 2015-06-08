import igraph
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import DataFrame, Series, read_csv

wd = '/Users/emanuel/Projects/projects/liverx'

# Import data-set
dataset = read_csv(wd + '/data/result_swath_v2.3.7_fRfS.csv').dropna()
dataset[dataset['adj.pvalue'] == 0] = dataset[dataset['adj.pvalue'] != 0]['adj.pvalue'].min()
dataset['log10.adj.pvalue'] = -np.log10(dataset['adj.pvalue'])
dataset_proteins = set(dataset['Protein'])

# Import PPI network
network = read_csv(wd + '/files/genemania_mouse_network_filtered.txt', sep='\t')
network = network[[a in dataset_proteins and b in dataset_proteins for a, b in zip(*(network['Gene_A'], network['Gene_B']))]]
network_nodes = set(network['Gene_A']).intersection(network['Gene_B'])
print '[INFO] Network imported: ', network.shape

# Condition
condition = 'B6_FASTED_T2-T0 - S9_FASTED_T2-T0'

sub_dataset = dataset[dataset['Label'] == condition].set_index('Protein')
sub_dataset = sub_dataset['log10.adj.pvalue'].to_dict()
sub_dataset = {k: sub_dataset[k] for k in sub_dataset if k in network_nodes}

# Export node file
nodes_files = Series(sub_dataset)
nodes_files.name = 'score1'
nodes_files.to_csv('/Users/emanuel/Projects/projects/liverx/src/liverx/nodes_weights.txt', header=True, index_label='#node', sep='\t')
print '[INFO] Node weights exported: ', nodes_files.shape

# Export edges file
sub_network = network[[a in sub_dataset or b in sub_dataset for a, b in zip(*(network['Gene_A'], network['Gene_B']))]]
sub_network['score1'] = 0
sub_network.columns = ['#nodeA', 'nodeB', 'score1']
sub_network.to_csv('/Users/emanuel/Projects/projects/liverx/src/liverx/edges_weights.txt', sep='\t', index=False)
print '[INFO] Edges weights exported: ', sub_network.shape