import subprocess
from liverx import wd

results_folder, fdr_thresholds = '%s/files/network_enrichment' % wd, [1e-3]

# ---- Iterate over the different FDR thresholds
for fdr_thres in fdr_thresholds:
    print '[INFO] Running BioNet with FDR threshold: %.0e' % fdr_thres

    # Run BioNet active modules
    subprocess.call('Rscript %s/src/liverx/network_active_modules.r %e %s/' % (wd, fdr_thres, results_folder), shell=True)