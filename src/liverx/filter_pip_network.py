import igraph
from pandas import read_csv, DataFrame
from liverx import wd

# ---- Create ID maps
# Id map between Ensemble and UniProt
id_map = read_csv(wd + '/files/10090_reviewed_uniprot_2_string.04_2015.tsv', sep='\t')
id_map = id_map.drop_duplicates(subset=['string_id'])
id_map = id_map.drop_duplicates(subset=['uniprot_ac|uniprot_id'])
id_map = id_map.set_index('string_id')['uniprot_ac|uniprot_id']
id_map.to_csv('%s/files/String_to_UniProt.tab' % wd, sep='\t')
id_map = id_map.to_dict()
print '[INFO] String to Uniprot map: ', len(id_map)

# ---- Import STRING mouse network
# Import network
info = read_csv(wd + '/files/10090.protein.links.v10.txt', sep=' ')
print '[INFO] Mouse String PPI network imported', len(info)

# Remove interactions with weight lower than X threshold
thres = 800

info = info[info['combined_score'] >= thres]
print '[INFO] Weight threshold %d: %d' % (thres, len(info))

# Remove String version from protein IDs
info['p1'] = [p.split('.')[1] for p in info['protein1']]
info['p2'] = [p.split('.')[1] for p in info['protein2']]

# Build igraph network
network = igraph.Graph(directed=False)
network.add_vertices(list(set(info['p1']).union(info['p2'])))
network.add_edges(list(zip(*(info['p1'], info['p2']))))
print '[INFO] Mouse PIP network', network.summary()

# Remove vertices without mapping ID
network = network.subgraph(set(id_map).intersection(network.vs['name']))
print '[INFO] Remove vertices without mapping ID: ', network.summary()

# Remove self-loops and multiple edges
network = network.simplify(True, True, 'first')
print '[INFO] Remove self-loops and multiple edges: ', network.summary()

# Create vertices uniprot id attribute
network.vs['uniprot'] = [id_map[v['name']].split('|')[1] for v in network.vs]

# Export to file
subnetwork = DataFrame([network.vs[(e.source, e.target)]['uniprot'] for e in network.es], columns=['protein1', 'protein2'])
subnetwork.to_csv('%s/files/string_mouse_network_filtered_%d.txt' % (wd, thres), sep='\t', index=False)
print '[INFO] Filtered network exported to: ', '%s/files/string_mouse_network_filtered_%d.txt' % (wd, thres)