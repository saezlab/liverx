from pandas import read_csv
from liverx import wd

# ---- Create ID maps
# Id map between Ensemble and UniProt
id_map = read_csv(wd + '/files/10090_reviewed_uniprot_2_string.04_2015.tsv', sep='\t')
id_map = id_map.drop_duplicates(subset=['string_id'])
id_map = id_map.drop_duplicates(subset=['uniprot_ac|uniprot_id'])
id_map = id_map.set_index('string_id')['uniprot_ac|uniprot_id']
id_map.to_csv('%s/files/String_to_UniProt.tab', sep='\t')
id_map = id_map.to_dict()
print '[INFO] String to Uniprot map: ', len(id_map)

# ---- Import STRING mouse network
# Import network
network = read_csv(wd + '/files/10090.protein.links.v10.txt', sep=' ')
print '[INFO] Mouse PIP network', network.shape

# Remove interactions with weight lower than X threshold
thres = 800
print '[INFO] Weight threshold: %d' % thres

subnetwork = network[network['combined_score'] >= thres]
print '[INFO] Remove interactions with weight below %d: ' % thres, subnetwork.shape

# Remove String version from protein IDs
subnetwork['protein1'] = [p.split('.')[1] for p in subnetwork['protein1']]
subnetwork['protein2'] = [p.split('.')[1] for p in subnetwork['protein2']]
print '[INFO] Protein IDs parsed: ', subnetwork.shape

# Remove interactions without mapping ID
subnetwork = subnetwork[[a in id_map and b in id_map for a, b in zip(*(subnetwork['protein1'], subnetwork['protein2']))]]
print '[INFO] Remove interactions without any measured node: ', subnetwork.shape

# Convert ENS to Uniprot
subnetwork['protein1'] = [id_map[i].split('|')[1] for i in subnetwork['protein1']]
subnetwork['protein2'] = [id_map[i].split('|')[1] for i in subnetwork['protein2']]

# Remove self interactions
subnetwork = subnetwork[[a != b for a, b in zip(*(subnetwork['protein1'], subnetwork['protein2']))]]
print '[INFO] Remove self interactions: ', subnetwork.shape

# Export to file
subnetwork[['protein1', 'protein2']].to_csv('%s/files/string_mouse_network_filtered_%d.txt' % (wd, thres), sep='\t', index=False)