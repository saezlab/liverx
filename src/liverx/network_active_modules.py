import pydot
import subprocess
from liverx import wd
from pandas import read_csv

results_folder, fdr_thresholds = '%s/files/network_enrichment' % wd, [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]

# ---- Iterate over the different FDR thresholds
for fdr_thres in fdr_thresholds:
    print '[INFO] Running BioNet with FDR threshold: %.0e' % fdr_thres

    # Run BioNet active modules
    subprocess.call('Rscript %s/src/liverx/network_active_modules.r %e %s/' % (wd, fdr_thres, results_folder), shell=True)

    for hypothesis in ['H2', 'H4']:
        # Import BioNet result network [hypothesis]_[FDR]_network.sif, e.g. H2_1e-07_network.xgmml
        graph = pydot.Dot(graph_type='digraph', rankdir='LR')
        graph_df = read_csv('%s/%s_%.0e_network.sif' % (results_folder, hypothesis, fdr_thres), sep='\t', header=None)

        for s, t in zip(*(graph_df[0], graph_df[2])):
            graph.add_edge(pydot.Edge(s, t))

        graph.write_pdf('%s/%s_%.0e_network.pdf' % (results_folder, hypothesis, fdr_thres))
        print '[INFO] Network exported: %s, %.0e' % (hypothesis, fdr_thres)