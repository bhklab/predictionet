### function evaluating the stability of the network given the data and method
## data: matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## categories: list of categories for each of the nodes (genes) in 'data' matrix. Categories should be integers from 1 to n. If method='bayesnet',  this parameter shoudl be specified by the user.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## priors: matrix of prior information avalable for gene-gene interaction (parents in columns, children in rows). Values may be probabilities or any other weights (citations count for instance)
## priors.count: 'TRUE' if priors specified by the user are number of citations for each interaction, 'FALSE' if probabilities are reported instead
## priors.weight: real value in [0, 1] specifying the weight to put on the priors (0=only the data are used, 1=only the priors are used to infer the topology of the network). If 'priors.weight' is missing it will be optimized gene by hene in an automatic way.
## maxparents: maximum number of parents allowed for each gene
## subset: vector of indices to select only subset of the observations
## method: "bayesnet" for bayesian network inference with the catnet package, "regrnet" for regression-based network inference
`net.stab` <- 
function(data, categories, perturbations, priors, priors.count=TRUE, priors.weight=0.5, maxparents=3, subset, method=c("bayesnet", "regrnet")) {
### evaluate stability by checking existence/absence of edges and direction of edges 
	
}
