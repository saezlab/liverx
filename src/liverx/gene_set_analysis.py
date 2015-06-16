import pickle
import os.path
import numpy as np
from statsmodels.stats.multitest import multipletests
from liverx import wd
from scipy.stats.distributions import hypergeom
from bioservices import KEGG, KEGGParser, QuickGO
from pandas import DataFrame, Series, read_csv

# Import ID maps
id_map_genename = read_csv('%s/files/GeneMania_Ensemble_GeneName.tab' % wd, sep='\t', header=None, index_col=1).to_dict()[0]
id_map_uniprot = read_csv('%s/files/GeneMania_Ensemble_UniProt.tab' % wd, sep='\t', header=None, index_col=1).to_dict()[0]

# Import network
network = read_csv('%s/files/genemania_mouse_network_filtered.txt' % wd, sep='\t')
network_genes = set(network['Gene_A']).intersection(network['Gene_B'])
network_genes_ens = {id_map_uniprot[i] for i in network_genes}

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

print '[INFO] GO Terms UniProt list dictionary built'

# ---- Import sub-network
hypothesis, fdr_thres = 'H4', 1e-07
sub_network = read_csv('%s/files/network_enrichment/test.sif' % wd, header=None, sep='\t')
sub_network_genes = {i for i in set(sub_network[0]).intersection(sub_network[2])}
sub_network_genes_ens = {id_map_uniprot[i] for i in sub_network_genes}

# ---- KEGG pathways enrichment analysis
# hypergeom.sf(x, M, n, N, loc=0)
# M: total number of objects,
# n: total number of type I objects
# N: total number of type I objects drawn without replacement
kegg_pathways_hyper = {p: hypergeom.sf(
    len(sub_network_genes_ens.intersection(kegg_pathways_genes[p])),
    len(network_genes_ens),
    len(network_genes_ens.intersection(kegg_pathways_genes[p])),
    len(sub_network_genes_ens)
) for p in kegg_pathways_genes if len(sub_network_genes_ens.intersection(kegg_pathways_genes[p])) > 0}
print '[INFO] KEGG pathways enrichment done'

# ---- GO terms enrichment analysis
# hypergeom.sf(x, M, n, N, loc=0)
# M: total number of objects,
# n: total number of type I objects
# N: total number of type I objects drawn without replacement
go_terms_hyper = {go: hypergeom.sf(
    len(sub_network_genes.intersection(go_terms[go])),
    len(network_genes),
    len(network_genes.intersection(go_terms[go])),
    len(sub_network_genes)
) for go in go_terms if len(sub_network_genes.intersection(go_terms[go])) > 0}
print '[INFO] GO terms enrichment done'