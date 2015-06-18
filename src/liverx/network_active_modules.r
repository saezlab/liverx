library(BioNet)
library(reshape)

# Configure workplace
setwd('~/Projects/projects/liverx/')

results_folder <- 'files/network_enrichment'
message(paste('Results folder set to:', results_folder))

# Import data-set
dataset <- read.csv('data/result_swath_v2.3.7_fRfS.csv')

# Import network
network <- loadNetwork.tab('files/genemania_mouse_network_filtered_1e-04.txt', header=T, directed=F, format='graphNEL')

# Import hypothesis conditions
comparisons <- read.table('files/hypothesis_conditions.txt', sep='\t', header=T)

# Sub-network for condition
lapply(c(1e-2, 1e-3, 1e-4), function (fdr_thres) {
  lapply(c('H2', 'H4'), function (condition) {
    # Get condition comparisons
    conditions <- comparisons[which(comparisons['hypothesis'] == condition), 'label']
    
    # Sub-set data-set comparisons to hypothesis
    sub_dataset <- dataset[which(dataset[,'Label'] %in% conditions),]
    sub_dataset <- na.omit(cast(sub_dataset[, c('Protein', 'Label', 'adj.pvalue')], Protein ~ Label))
    
    proteins <- sub_dataset[,'Protein']
    sub_dataset <-  as.matrix(sub_dataset[, -1])
    rownames(sub_dataset) <- proteins
    
    pvals <- aggrPvals(sub_dataset, order=dim(sub_dataset)[2], plot=F)
    pvals[pvals == 0] <- min(pvals[pvals != 0])
    
    # Reduce PPI network size
    sub_network <- subNetwork(names(pvals), network)
    sub_network <- largestComp(sub_network)
    sub_network <- rmSelfLoops(sub_network)
    
    # Fit Beta-uniform mixture model
    fb <- fitBumModel(pvals, plot=F)
    scores <- scoreNodes(network=sub_network, fb=fb, fdr=fdr_thres)
    scores <- round(scores, 4)
    
    # Write files to run Heinz
    edges_file <- paste(results_folder, condition, '_', fdr_thres, '_edges.txt', sep='')
    nodes_file <- paste(results_folder, condition, '_', fdr_thres, '_nodes.txt', sep='')
    result_file <- paste(results_folder, condition, '_', fdr_thres, '_solution.txt', sep='')
    network_file <- paste(results_folder, condition, '_', fdr_thres, '_network', sep='')
    
    writeHeinzEdges(network=sub_network, file=edges_file, use.score=FALSE)
    writeHeinzNodes(network=sub_network, file=nodes_file, node.scores=scores)
    
    # Execute Heinz
    run_heinz_cmd <- paste('heinz -e ', edges_file, ' -n ', nodes_file, ' -o ', result_file, ' -t 3600',sep='')
    system(run_heinz_cmd)
    
    # Importing the solution of a maxinum-scoring subnetwork
    module <- readHeinzGraph(node.file=result_file, network=sub_network)
    
    nodeDataDefaults(module, attr='score') <- ''
    nodeData(module, n=nodes(module), attr='score') <- scores[nodes(module)]
    
    # Export network
    saveNetwork(module, file=network_file, type='sif')
  })
})