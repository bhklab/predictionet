### function fitting a bayesian network model
## data: Matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## priors: matrix of prior information avalable for gene-gene interaction (parents in columns, children in rows). Values may be probabilities or any other weights (citations count for instance)
## maxparents: maximum number of parents allowed for each gene
## seed: seed to get deterministic results
### returns a bayesian network model
'.fit.catnet' <- 
function(data, categories, perturbations, priors, priors.weight, maxparents=3, bayesnet.maxcomplexity=0, bayesnet.maxiter=100, seed) {
	#require(catnet)
	## perturbations should be a matrix of [0, 1] for the catnet package
	pp <- matrix(0, nrow=nrow(perturbations), ncol=ncol(perturbations), dimnames=dimnames(perturbations))
	pp[perturbations] <- 1
	pertwurbations <- pp
	if(!missing(seed)) { catnet::cnSetSeed(seed) }
	## be aware that catnet package consider any adjacency matrix to have parents in COLUMNS and children in ROWS, that is the inverse of the predictionet package
	#edgrel <- matrix(0, nrow=ncol(data), ncol=ncol(data), dimnames=list(colnames(data), colnames(data)))
	## topology identified from the priors
	priorparents <- sapply(1:ncol(priors), function(x, y) { tt <- which(y[ , x] > 0.5); if(length(tt) < 1) { tt <- NULL }; return(list(tt)); }, y=adj.remove.cycles(adjmat=priors != 0.5)[[1]])
	names(priorparents) <- colnames(priors)
	if(priors.weight == 1) {
		## use only priors to fix the topology
		## check if the priors are acyclic
		## create network
		ee <- catnet::cnNew(nodes=colnames(data), cats=categories, parents=priorparents)
		#ee <- catnet::cnSetProb(object=ee, data=t(data), nodeCats=categories, perturbations=t(perturbations))
		ee <- catnet::cnSetProb(object=ee, data=t(data), perturbations=t(perturbations))
		#ee <- .cnUpdateCat(object=ee, cats=categories)
		myparents <- lapply(ee@parents, function(x, y) { if(!is.null(x)) { x <- y[x]}; return(x); }, y=ee@nodes)
		names(myparents) <- ee@nodes
		return(list("varnames"=colnames(data), "input"=myparents, "model"=ee, "edge.relevance"=edgerel))
	} else {
		## use data and priors infer the best topology
		if(priors.weight == 0) {
			## do not use the priors for identifying the network topology
			## uniform priors
			priors[] <- 0.5
		}
		## seed the network inference with the prior topology
		priororder <- catnet::cnOrder(priorparents)
		ee.prior <- catnet::cnSearchOrder(data=t(data), perturbations=t(perturbations), maxParentSet=maxparents, nodeOrder=priororder, edgeProb=t(priors), maxComplexity=bayesnet.maxcomplexity)
		#ee.prior <- NULL
		ee <- catnet::cnSearchSA(data=t(data), nodeCats=categories, perturbations=t(perturbations), selectMode="BIC", maxParentSet=maxparents, priorSearch=ee.prior, edgeProb=t(priors), echo=FALSE, maxComplexity=bayesnet.maxcomplexity, maxIter=bayesnet.maxiter)
		ee <- catnet::cnFindBIC(ee)
		myparents <- lapply(ee@parents, function(x, y) { if(!is.null(x)) { x <- y[x]}; return(x); }, y=ee@nodes)
		names(myparents) <- ee@nodes
		## compute edge relevance as the likelihood of the target/child node
		#tt <- exp(sapply(1:ncol(data), function(x, y, z, w) { return(catnet::cnNodeLoglik(object=y, node=x, data=z, perturbations=w)) }, y=ee, z=t(data), w=t(perturbations)))
		#names(tt) <- colnames(data)
		#edgerel <- mapply(function(x,y,z) { rr <- rep(0, length(z)); names(rr) <- z; rr[x] <- y; return(rr) }, x=ee@parents, y=tt, MoreArgs=list(z=colnames(data)))
		#dimnames(edgerel) <- list(colnames(data), colnames(data))
		edgerel <- matrix(NA, nrow=ncol(data), ncol=ncol(data), dimnames=list(colnames(data), colnames(data)))
		
		return(list("varnames"=colnames(data), "input"=myparents, "model"=ee, "edge.relevance"=edgerel))
	}
}
