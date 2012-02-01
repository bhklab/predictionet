### function selecting at most maxparents parents for each target variable
## data: Matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## priors: matrix of prior information avalable for gene-gene interaction (parents in columns, children in rows). Values may be probabilities or any other weights (citations count for instance)
## predn: indices or names of genes to fit during network inference. If missing, all the genes will be used for network inference
## priors.weight: real value in [0, 1] specifying the weight to put on the priors (0=only the data are used, 1=only the priors are used to infer the topology of the network). If 'priors.weight' is missing it will be optimized gene by hene in an automatic way.
## maxparents: maximum number of parents allowed for each gene
### returns a list of local networks (myparents), the causality score matrix (edge.relevance) and the names for all variables in the network (varnames)
`.fit.regrnet.causal2` <- 
function(data, perturbations, priors, predn, maxparents=3, priors.weight=0.5, causal=TRUE, seed) {
	if(causal && maxparents > (ncol(data) * 0.5)) { warning("maximum number of parents may be too large, causal inference requires sparsity in the inferred network; please decrease maxparents parameter for better results!") }
	if(!missing(seed)) { set.seed(seed) }
	if(missing(predn) || is.null(predn) || length(predn) == 0) { predn <- 1:ncol(data) } else { if(is.character(predn)) { predn <- match(predn, dimnames(data)[[2]]) } else { if(!is.numeric(predn) || !all(predn %in% 1:ncol(data))) { stop("parameter 'predn' should contain either the names or the indices of the variables to predict!")} } }
	nn <- dimnames(data)[[2]][predn] ## variables (genes) to fit during network inference
	regrnet <- NULL
	
	########################
	### data preparation
	########################
	dd <- data.frame(data)
	diag(priors) <- 0 ## no self-loops
	
	########################
	### find ranked parents for all genes
	########################
	rrnetw <- .rank.genes.causal.perturbations2(priors=priors, data=dd, perturbations=perturbations, predn=predn, priors.weight=priors.weight, maxparents=maxparents, causal=causal)
	mat.ranked.parents <- rrnetw$parents
	myparents <- lapply(apply(mat.ranked.parents, 1, function(x) { x <- x[!is.na(x)]; if(length(x) == 0) { x <- NULL };  return(list(x)); }), function(x) { return(x[[1]]) })
	

	return(list("varnames"=colnames(data), "input"=myparents, "edge.relevance"=rrnetw$score))
}
