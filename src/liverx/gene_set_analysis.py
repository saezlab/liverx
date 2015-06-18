import re
import pydot
import pickle
import os.path
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import rcParams
from liverx import wd
from statsmodels.stats.multitest import multipletests
from scipy.stats.distributions import hypergeom
from bioservices import KEGG, KEGGParser, QuickGO
from pandas import DataFrame, Series, read_csv

sns.set_style('white')

# Import ID maps
id_map_genename = read_csv('%s/files/GeneMania_Ensemble_GeneName.tab' % wd, sep='\t', header=None, index_col=1).to_dict()[0]
id_map_uniprot = read_csv('%s/files/GeneMania_Ensemble_UniProt.tab' % wd, sep='\t', header=None, index_col=1).to_dict()[0]

# Import network
network = read_csv('%s/files/genemania_mouse_network_filtered.txt' % wd, sep='\t')
network_genes = set(network['Gene_A']).intersection(network['Gene_B'])
network_genes_ens = {id_map_uniprot[i] for i in network_genes}

# ---- Set-up QuickGO bioservice
quickgo = QuickGO(cache=True)

# ---- Set-up KEGG bioservice
kegg, kegg_parser = KEGG(cache=True), KEGGParser()

kegg.organism = 'mmu'
print '[INFO] KEGG service configured'

kegg_pathways = {p: kegg.parse_kgml_pathway(p) for p in kegg.pathwayIds}
print '[INFO] KEGG pathways extracted: ', len(kegg_pathways)

# Convert KEGG pathways Gene Name to Ensemble
kegg_pathways_genes = {p: {i['gene_names'].split(', ')[0] for i in kegg_pathways[p]['entries'] if i['type'] == 'gene'} for p in kegg_pathways}
kegg_pathways_genes = {p: {id_map_genename[e] for e in kegg_pathways_genes[p] if e in id_map_genename} for p in kegg_pathways_genes}

# ---- Set-up GO Terms gene list
go_terms_file = '%s/files/go_terms_uniprot.pickle' % wd

if os.path.isfile(go_terms_file):
    with open(go_terms_file, 'rb') as handle:
        go_terms = pickle.load(handle)

else:
    go_terms = read_csv('%s/files/gene_association.goa_ref_mouse' % wd, sep='\t', skiprows=12, header=None)[[4, 10]]
    go_terms_set = set(go_terms[4])
    go_terms = {go: set(go_terms.loc[go_terms[4] == go, 10]) for go in go_terms_set}
    go_terms = {go: {p.split('|')[0] for p in go_terms[go]} for go in go_terms}

    with open(go_terms_file, 'wb') as handle:
        pickle.dump(go_terms, handle)

print '[INFO] GO Terms UniProt list dictionary built: ', len(go_terms)

go_terms = {go: go_terms[go] for go in go_terms if '<namespace>biological_process</namespace>' in quickgo.Term(go)}
print '[INFO] Consider only biological process GO Terms: ', len(go_terms)

for hypothesis in ['H2', 'H4']:
    # ---- Import sub-network
    fdr_thres = '0.01'
    sub_network = read_csv('%s/files/network_enrichment/%s_%s_network.sif' % (wd, hypothesis, fdr_thres), header=None, sep='\t')
    sub_network_genes = {i for i in set(sub_network[0]).union(sub_network[2])}
    sub_network_genes_ens = {id_map_uniprot[i] for i in sub_network_genes}

    # ---- KEGG pathways enrichment analysis
    # hypergeom.sf(x, M, n, N, loc=0)
    # M: total number of objects,
    # n: total number of type I objects
    # N: total number of type I objects drawn without replacement
    kegg_pathways_hyper = {p: (hypergeom.sf(
        len(sub_network_genes_ens.intersection(kegg_pathways_genes[p])),
        len(network_genes_ens),
        len(network_genes_ens.intersection(kegg_pathways_genes[p])),
        len(sub_network_genes_ens)
    ), len(sub_network_genes_ens.intersection(kegg_pathways_genes[p]))) for p in kegg_pathways_genes if len(sub_network_genes_ens.intersection(kegg_pathways_genes[p])) > 0}
    print '[INFO] KEGG pathways enrichment done'

    kegg_pathways_hyper = DataFrame(kegg_pathways_hyper, index=['pvalue', 'intersection']).T
    kegg_pathways_hyper['adj.pvalue'] = multipletests(kegg_pathways_hyper['pvalue'], method='fdr_bh')[1]
    kegg_pathways_hyper = kegg_pathways_hyper.sort('adj.pvalue', ascending=False)
    kegg_pathways_hyper = kegg_pathways_hyper[kegg_pathways_hyper['adj.pvalue'] < 0.05]
    kegg_pathways_hyper['name'] = [re.findall('NAME\s*(.*) - Mus musculus\n?', kegg.get(p))[0] for p in kegg_pathways_hyper.index]
    kegg_pathways_hyper = kegg_pathways_hyper[kegg_pathways_hyper['adj.pvalue'] != 0.0]

    rcParams['figure.figsize'] = 5, 11
    y_pos = np.arange(len(kegg_pathways_hyper))
    plt.barh(y_pos, -np.log10(kegg_pathways_hyper['adj.pvalue']), lw=0, alpha=.8, align='center')
    sns.despine()
    plt.yticks(y_pos, kegg_pathways_hyper['name'])
    plt.xlabel('p-value (FDR -log10)')
    plt.title('KEGG mouse pathways (hyper-geometric test)')
    plt.savefig('%s/reports/Figure7_%s_kegg_enrichment.pdf' % (wd, hypothesis), bbox_inches='tight')
    plt.close('all')
    print '[INFO] KEGG pathways enrichment plotted'

    # ---- GO terms enrichment analysis
    # hypergeom.sf(x, M, n, N, loc=0)
    # M: total number of objects,
    # n: total number of type I objects
    # N: total number of type I objects drawn without replacement
    go_terms_hyper = {go: (hypergeom.sf(
        len(sub_network_genes.intersection(go_terms[go])),
        len(network_genes),
        len(network_genes.intersection(go_terms[go])),
        len(sub_network_genes)
    ), len(sub_network_genes.intersection(go_terms[go]))) for go in go_terms if len(sub_network_genes.intersection(go_terms[go])) > 0}
    print '[INFO] GO terms enrichment done'

    go_terms_hyper = DataFrame(go_terms_hyper, index=['pvalue', 'intersection']).T.dropna()
    go_terms_hyper['adj.pvalue'] = multipletests(go_terms_hyper['pvalue'], method='fdr_bh')[1]
    go_terms_hyper = go_terms_hyper.sort('adj.pvalue', ascending=False)
    go_terms_hyper = go_terms_hyper[go_terms_hyper['adj.pvalue'] < 0.05]
    go_terms_hyper['name'] = [re.findall('name: (.*)\n?', quickgo.Term(go, frmt='obo'))[0] for go in go_terms_hyper.index]
    go_terms_hyper = go_terms_hyper[go_terms_hyper['intersection'] > 2]
    go_terms_hyper = go_terms_hyper[go_terms_hyper['adj.pvalue'] != 0.0]

    rcParams['figure.figsize'] = 5, 0.25 * len(go_terms_hyper)
    y_pos = np.arange(len(go_terms_hyper))
    plt.barh(y_pos, -np.log10(go_terms_hyper['adj.pvalue']), lw=0, alpha=.8, align='center')
    sns.despine()
    plt.yticks(y_pos, go_terms_hyper['name'])
    plt.xlabel('p-value (FDR -log10)')
    plt.title('GO Terms (hyper-geometric test)')
    plt.savefig('%s/reports/Figure7_%s_goterms_enrichment.pdf' % (wd, hypothesis), bbox_inches='tight')
    plt.close('all')
    print '[INFO] GO terms enrichment plotted'

    # ---- Plot sub-network
    graph = pydot.Dot(graph_type='graph', rankdir='LR')

    for s, t in zip(*(sub_network[0], sub_network[2])):
        graph.add_edge(pydot.Edge(s.split('_')[0], t.split('_')[0]))

    graph.write_pdf('%s/reports/Figure7_%s_active_subnetwork.pdf' % (wd, hypothesis))
    print '[INFO] Network exported: %s, %s' % (hypothesis, fdr_thres)