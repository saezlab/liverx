import re
import pickle
import igraph
import os.path
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import rcParams
from liverx import wd
from matplotlib.colors import rgb2hex
from statsmodels.stats.multitest import multipletests
from scipy.stats.distributions import hypergeom
from bioservices import KEGG, KEGGParser, QuickGO, UniProt
from pandas import DataFrame, read_csv

sns.set(style='ticks', palette='pastel', color_codes=True)

# ---- Import network
network = read_csv('%s/files/string_mouse_network_filtered_800.txt' % wd, sep='\t')
network_proteins = set(network['protein1']).intersection(network['protein2'])

# ---- Set-up UniProt
uniprot = UniProt(cache=True)

# ---- Set-up QuickGO bioservice
quickgo = QuickGO(cache=True)

# ---- Set-up KEGG bioservice
kegg, kegg_parser = KEGG(cache=True), KEGGParser()

kegg.organism = 'mmu'
print '[INFO] KEGG service configured'

kegg_pathways = {p: kegg.parse_kgml_pathway(p) for p in kegg.pathwayIds}
print '[INFO] KEGG pathways extracted: ', len(kegg_pathways)

# Convert KEGG pathways Gene Name to UniProt
k2u = kegg.conv('uniprot', 'mmu')

kegg_pathways_proteins = {p: {k2u[x].split(':')[1] for i in kegg_pathways[p]['entries'] if i['type'] == 'gene' for x in i['name'].split(' ') if x in k2u} for p in kegg_pathways}

kegg_uniprot_acc_map = {x for p in kegg_pathways_proteins for x in kegg_pathways_proteins[p]}
kegg_uniprot_acc_map = {p: uniprot.get_fasta(str(p)).split(' ')[0].split('|')[2] for p in kegg_uniprot_acc_map}

kegg_pathways_proteins = {p: {kegg_uniprot_acc_map[i] for i in kegg_pathways_proteins[p]} for p in kegg_pathways_proteins}
print '[INFO] KEGG pathways Ids converted to UniProt: ', len(kegg_pathways_proteins)

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

# ---- Run Hypergeometric test
for hypothesis, fdr_thres in [('H2', '0.05'), ('H4', '0.05')]:

    # ---- Import sub-network
    subnetwork = read_csv('%s/files/network_enrichment/%s_%s_network.sif' % (wd, hypothesis, fdr_thres), sep='\t', names=['p1', 'i', 'p2'])

    subnetwork_i = igraph.Graph(directed=False)
    subnetwork_i.add_vertices(list(set(subnetwork['p1']).union(subnetwork['p2'])))
    subnetwork_i.add_edges([(p1, p2) for p1, p2 in subnetwork[['p1', 'p2']].values])
    print '[INFO] String network: ', subnetwork_i.summary()

    # Generate network layout positions
    palette = [rgb2hex((r, g, b)) for r, g, b in sns.color_palette('Paired')[:2]]
    layout = subnetwork_i.layout_fruchterman_reingold(maxiter=10000, area=50 * (len(subnetwork_i.vs) ** 2))
    print '[INFO] String network layout created: ', subnetwork_i.summary()

    subnetwork_proteins = set(subnetwork_i.vs['name'])

    # ---- KEGG pathways enrichment analysis
    # hypergeom.sf(x, M, n, N, loc=0)
    # M: total number of objects,
    # n: total number of type I objects
    # N: total number of type I objects drawn without replacement
    kegg_pathways_hyper = {p: (hypergeom.sf(
        len(subnetwork_proteins.intersection(kegg_pathways_proteins[p])),
        len(network_proteins),
        len(network_proteins.intersection(kegg_pathways_proteins[p])),
        len(subnetwork_proteins)
    ), len(subnetwork_proteins.intersection(kegg_pathways_proteins[p]))) for p in kegg_pathways_proteins if len(subnetwork_proteins.intersection(kegg_pathways_proteins[p])) > 0}
    print '[INFO] KEGG pathways enrichment done'

    kegg_pathways_hyper = DataFrame(kegg_pathways_hyper, index=['pvalue', 'intersection']).T
    kegg_pathways_hyper['adj.pvalue'] = multipletests(kegg_pathways_hyper['pvalue'], method='fdr_bh')[1]
    kegg_pathways_hyper = kegg_pathways_hyper.sort('adj.pvalue', ascending=False)
    kegg_pathways_hyper = kegg_pathways_hyper[kegg_pathways_hyper['adj.pvalue'] < 0.05]
    kegg_pathways_hyper['name'] = [re.findall('NAME\s*(.*) - Mus musculus\n?', kegg.get(p))[0] for p in kegg_pathways_hyper.index]
    kegg_pathways_hyper = kegg_pathways_hyper[kegg_pathways_hyper['adj.pvalue'] != 0.0]

    # Plot PPI network of the Kegg pathways
    for set_id, set_name in kegg_pathways_hyper['name'].tail(10).to_dict().items():
        subnetwork_i.vs['label'] = [n.split('_')[0] if n in kegg_pathways_proteins[set_id] else '' for n in subnetwork_i.vs['name']]
        subnetwork_i.vs['shape'] = ['circle' for n in subnetwork_i.vs['name']]
        subnetwork_i.vs['color'] = [palette[1] if n in kegg_pathways_proteins[set_id] else palette[0] for n in subnetwork_i.vs['name']]
        subnetwork_i.vs['size'] = [17 if n in kegg_pathways_proteins[set_id] else 7 for n in subnetwork_i.vs['name']]

        igraph.plot(
            subnetwork_i,
            layout=layout,
            vertex_label_size=5,
            vertex_label_color='white',
            bbox=(0, 0, 360, 360),
            vertex_frame_width=0,
            edge_width=.5,
            edge_color='#ececec',
            target='%s/reports/Figure7_%s_%s_kegg_enrichment_network_%s_%s.pdf' % (wd, hypothesis, fdr_thres, set_id, set_name)
        )

    # Plot Kegg pathways enrichment values barplot
    rcParams['figure.figsize'] = 5, 0.3 * len(kegg_pathways_hyper)
    colours, y_pos = sns.color_palette('Paired', 2), [x + 1.5 for x in range(len(kegg_pathways_hyper))]

    plt.barh(y_pos, -np.log10(kegg_pathways_hyper['pvalue']), lw=0, align='center', height=.5, color=colours[0], label='p-value')
    plt.barh(y_pos, -np.log10(kegg_pathways_hyper['adj.pvalue']), lw=0, align='center', height=.5, color=colours[1], label='FDR')
    plt.yticks(y_pos, kegg_pathways_hyper['name'])

    plt.axvline(-np.log10(0.05), ls='--', lw=0.4, c='gray')
    plt.axvline(-np.log10(0.01), ls='--', lw=0.4, c='gray')

    plt.text(-np.log10(0.05) * 1.01, .5, '5%', ha='left', color='gray', fontsize=9)
    plt.text(-np.log10(0.01) * 1.01, .5, '1%', ha='left', color='gray', fontsize=9)

    sns.despine()
    plt.xlabel('-log10')
    plt.title('KEGG pathways enrichment')
    plt.legend(loc=4)
    plt.savefig('%s/reports/Figure7_%s_%s_kegg_enrichment.pdf' % (wd, hypothesis, fdr_thres), bbox_inches='tight')
    plt.close('all')
    print '[INFO] KEGG pathways enrichment plotted'

    # ---- GO terms enrichment analysis
    # hypergeom.sf(x, M, n, N, loc=0)
    # M: total number of objects,
    # n: total number of type I objects
    # N: total number of type I objects drawn without replacement
    go_terms_hyper = {go: (hypergeom.sf(
        len(subnetwork_proteins.intersection(go_terms[go])),
        len(network_proteins),
        len(network_proteins.intersection(go_terms[go])),
        len(subnetwork_proteins)
    ), len(subnetwork_proteins.intersection(go_terms[go]))) for go in go_terms if len(subnetwork_proteins.intersection(go_terms[go])) > 0}
    print '[INFO] GO terms enrichment done'

    go_terms_hyper = DataFrame(go_terms_hyper, index=['pvalue', 'intersection']).T.dropna()
    go_terms_hyper['adj.pvalue'] = multipletests(go_terms_hyper['pvalue'], method='fdr_bh')[1]
    go_terms_hyper = go_terms_hyper.sort('adj.pvalue', ascending=False)
    go_terms_hyper = go_terms_hyper[go_terms_hyper['adj.pvalue'] < 0.05]
    go_terms_hyper['name'] = [re.findall('name: (.*)\n?', quickgo.Term(go, frmt='obo'))[0] for go in go_terms_hyper.index]
    go_terms_hyper = go_terms_hyper[go_terms_hyper['intersection'] > 2]
    go_terms_hyper = go_terms_hyper[go_terms_hyper['adj.pvalue'] != 0.0]

    # Plot PPI network of the GO terms
    for set_id, set_name in go_terms_hyper['name'].tail(10).to_dict().items():
        subnetwork_i.vs['label'] = [n.split('_')[0] if n in go_terms[set_id] else '' for n in subnetwork_i.vs['name']]
        subnetwork_i.vs['shape'] = ['circle' for n in subnetwork_i.vs['name']]
        subnetwork_i.vs['color'] = [palette[1] if n in go_terms[set_id] else palette[0] for n in subnetwork_i.vs['name']]
        subnetwork_i.vs['size'] = [17 if n in go_terms[set_id] else 7 for n in subnetwork_i.vs['name']]

        igraph.plot(
            subnetwork_i,
            layout=layout,
            vertex_label_size=5,
            vertex_label_color='white',
            bbox=(0, 0, 360, 360),
            vertex_frame_width=0,
            edge_width=.5,
            edge_color='#ECECEC',
            target='%s/reports/Figure7_%s_%s_goterms_enrichment_network_%s_%s.pdf' % (wd, hypothesis, fdr_thres, set_id, set_name)
        )

    # Plot GO terms enrichment values barplot
    rcParams['figure.figsize'] = 5, 0.3 * len(go_terms_hyper)
    colours, y_pos = sns.color_palette('Paired', 2), [x + 1.5 for x in range(len(go_terms_hyper))]

    plt.barh(y_pos, -np.log10(go_terms_hyper['pvalue']), lw=0, align='center', height=.5, color=colours[0], label='p-value')
    plt.barh(y_pos, -np.log10(go_terms_hyper['adj.pvalue']), lw=0, align='center', height=.5, color=colours[1], label='FDR')
    plt.yticks(y_pos, go_terms_hyper['name'])

    plt.axvline(-np.log10(0.05), ls='--', lw=0.4, c='gray')
    plt.axvline(-np.log10(0.01), ls='--', lw=0.4, c='gray')

    plt.text(-np.log10(0.05) * 1.01, .5, '5%', ha='left', color='gray', fontsize=9)
    plt.text(-np.log10(0.01) * 1.01, .5, '1%', ha='left', color='gray', fontsize=9)

    sns.despine()
    plt.xlabel('-log10')
    plt.title('KEGG pathways enrichment')
    plt.legend(loc=4)
    plt.savefig('%s/reports/Figure7_%s_%s_goterms_enrichment.pdf' % (wd, hypothesis, fdr_thres), bbox_inches='tight')
    plt.close('all')
    print '[INFO] GO terms enrichment plotted'
