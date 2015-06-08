from pandas import read_csv

wd = '/Users/emanuel/Projects/projects/liverx'

# Import identifiers converter
id_map = read_csv(wd + '/files/genemania_mouse_identifier_mappings.txt', sep='\t')
id_map = id_map[id_map['Source'] == 'Uniprot ID']
id_map = id_map[[x.endswith('_MOUSE') for x in id_map['Name']]]
id_map = id_map.drop_duplicates(subset=['Preferred_Name'])
id_map = id_map.set_index('Preferred_Name')['Name'].to_dict()
print '[INFO] ENS to Uniprot map: ', len(id_map)

# Import network
network = read_csv(wd + '/files/genemania_mouse_network.txt', sep='\t')
print '[INFO] GeneMania mouse PIP network', network.shape

# Remove interactions with weight lower than X threshold
thres = 10e-4
network = network[network['Weight'] > thres]
print '[INFO] Remove interactions with weight below %f: ' % thres, network.shape

# Remove self interactions
network = network[[a != b for a, b in zip(*(network['Gene_A'], network['Gene_B']))]]
print '[INFO] Remove self interactions: ', network.shape

# Remove interactions without any measured node
network = network[[a in id_map and b in id_map for a, b in zip(*(network['Gene_A'], network['Gene_B']))]]
print '[INFO] Remove interactions without any measured node: ', network.shape

# Convert ENS to Uniprot
network['Gene_A'] = [id_map[i] for i in network['Gene_A']]
network['Gene_B'] = [id_map[i] for i in network['Gene_B']]

# Export to file
network[['Gene_A', 'Gene_B']].to_csv(wd + '/files/genemania_mouse_network_filtered.txt', sep='\t', index=False)