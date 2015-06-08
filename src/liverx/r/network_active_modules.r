library(BioNet)
library(DLBCL)
library(reshape)

# Configure workplace
setwd('~/Projects/projects/liverx/src/liverx/r/')

# Import data-set
dataset <- read.csv('~/Projects/projects/liverx/data/result_swath_v2.3.7_fRfS.csv')

# Import network
network <- loadNetwork.tab('~/Projects/projects/liverx/files/genemania_mouse_network_filtered.txt', header=T, directed=F, format='graphNEL')

# Import hypothesis conditions
comparisons <- read.table('~/Projects/projects/liverx/files/hypothesis_conditions.txt', sep='\t', header=T)

# Sub-network for condition
lapply(c('H2', 'H4'), function (condition) {
  # Get condition comparisons
  conditions <- comparisons[which(comparisons['hypothesis'] == condition), 'label']
  
  # Sub-set data-set comparisons to hypothesis
  sub_dataset <- na.omit(dataset[which(dataset[,'Label'] %in% conditions),])
  sub_dataset <- sub_dataset[, c('Protein', 'Label', 'adj.pvalue')]
  sub_dataset <- na.omit(cast(sub_dataset, Protein ~ Label))
  rownames(sub_dataset) <- sub_dataset[,'Protein']
  sub_dataset <- sub_dataset[, -1]
  
  pvals <- aggrPvals(as.matrix(sub_dataset), order=length(sub_dataset), plot=F)
  names(pvals) <- rownames(sub_dataset)
  pvals[pvals == 0] <- min(pvals[pvals != 0])
  
  # Reduce PPI network size
  sub_network <- subNetwork(names(pvals), network)
  sub_network <- largestComp(sub_network)
  sub_network <- rmSelfLoops(sub_network)
  
  # Fit Beta-uniform mixture model
  fb <- fitBumModel(pvals)
  scores <- scoreNodes(network=sub_network, fb=fb, fdr=1e-8)
  plot(scores, pvals[names(scores)])
  
  # Write files to run Heinz
  writeHeinzEdges(network=sub_network, file=paste(condition, '_edges', sep=''), use.score=FALSE)
  writeHeinzNodes(network=sub_network, file=paste(condition, '_nodes', sep=''), node.scores=scores)
  
  # Execute Heinz
  run_heinz_cmd <- paste('./heinz -e ', condition, '_edges.txt -n ', condition, '_nodes.txt -o ', condition, '_solution.txt', sep='')
  system(run_heinz_cmd)
  
  # Importing the solution of a maxinum-scoring subnetwork
  module <- readHeinzGraph(node.file=paste(condition, '_solution.txt', sep=''), network=sub_network)
  
  nodeDataDefaults(module, attr='score') <- ''
  nodeData(module, n=nodes(module), attr='score') <- scores[nodes(module)]
  
  # Export network
  saveNetwork(module, file=paste(condition, '_module', sep=''), type='XGMML')
})