### Function inferring a network
## data: matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## categories: Either an integer specifying the number of categories used to discretize the variables or a list of categories for each of the variables in 'data' matrix. Categories should be integers from 1 to n. If method='bayesnet', this parameter should be specified by the user.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## priors: matrix of prior information avalable for gene-gene interaction (parents in columns, children in rows). Values may be probabilities or any other weights (citations count for instance)
## predn: indices or names of genes to fit during network inference. If missing, all the genes will be used for network inference
## priors.count: 'TRUE' if priors specified by the user are number of citations for each interaction, 'FALSE' if probabilities are reported instead
## priors.weight: real value in [0, 1] specifying the weight to put on the priors (0=only the data are used, 1=only the priors are used to infer the topology of the network). If 'priors.weight' is missing it will be optimized gene by hene in an automatic way.
## maxparents: maximum number of parents allowed for each gene
## subset: vector of indices to select only subset of the observations
## method: "bayesnet" for bayesian network inference with the catnet package, "regrnet" for regression-based network inference
## regrmodel: type of regression model to fit when 'regrnet' method is selected
## seed: seed to make the function fully deterministic
## bayesnet.maxcomplexity: maximum complexity for bayesian network inference
## bayesnet.maxiter: maximum number of iterations for bayesian network inference

`netinf` <-  function(data, categories, perturbations, priors, predn, priors.count=TRUE, priors.weight=0.5, maxparents=3, maxparents.push=FALSE, subset, method=c("regrnet", "bayesnet"),ensemble=FALSE, regrmodel=c("linear", "linear.penalized"), ensemble.model=c("full","best"), ensemble.maxnsol=3, causal=TRUE, seed=54321, bayesnet.maxcomplexity=0, bayesnet.maxiter=100) {
## select more genes to find a number of cauqal variables close or more than maxparents
	if(ncol(data) < 2) { stop("Number of variables is too small to infer a network!") }
	if(!missing(predn) && !is.null(predn) && (length(predn) < 2 && method!="regrnet.ensemble")) { stop("length of parameter 'predn' should be >= 2!") }
	maxparents <- ifelse(maxparents >= ncol(data), ncol(data)-1, maxparents)
	maxparents <- ifelse(maxparents < 2, 2, maxparents)
	method <- match.arg(method)
	regrmodel <- match.arg(regrmodel)
	ensemble.model <- match.arg(ensemble.model)
	if(missing(perturbations)) {
## create matrix of no perturbations
		perturbations <- matrix(0, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	}
	if(missing(priors)) {
## create a matrix of priors count = 0
		priors <- matrix(0, nrow=ncol(data), ncol=ncol(data), dimnames=list(dimnames(data)[[2]], dimnames(data)[[2]]))
	}
	if(priors.weight < 0 || priors.weight > 1) { stop("'priors.weight' should be in [0,1]!") }
	if(!missing(subset)) {
		data <- data[subset, , drop=FALSE]
		perturbations <- perturbations[subset, , drop=FALSE]
	}
	switch(method, 
		   "bayesnet"={
		   if(ensemble){
		   stop("ensemble bayesian network inference is not implemented yet!")
		   }else{
## fit bayesian network model
## priors
		   if(priors.count) {
## transform priors into probabilities
		   if(priors.weight == 0) { priors[] <- 0.5 } else {
## rescale the priors count so they lay in [0,1] with values more different from 0.5 propotionally to priors.weight
## maximum number of citations, without taking into the 5% largest counts
		   ma <- quantile(x=abs(priors[priors != 0]), probs=0.95)
		   pp <- (abs(priors) / ma)
		   pp[pp >= 1] <- 0.9
		   priors <- 0.5 + sign(priors) %*% (apply(X=pp, MARGIN=c(1, 2), FUN=function(x, y) { if(y < 0.5) { x <- ((1-y^x)) + x * (y) } else { x <- ((1-y^x)) + x * (1-y) }; return(x); }, y=1-priors.weight)) / 2
## FORMULA IS WRONG, NEED TO UPDATE IT
		   priors[priors >= 1] <- 0.99
		   }
		   } else { if(!all(priors >= 0 & priors <= 1)) { stop("if 'priors.count' is FALSE method 'bayesnet' requires priors to be specified as probabilities!") } }
## categories
		   if(!missing(categories)) {
## discretize gene expression data
		   if(is.numeric(categories)) {
		   if(length(categories) == 1) {
## use the same number of categories for each variable
		   categories <- rep(categories, ncol(data))
		   names(categories) <- dimnames(data)[[2]]
		   }
		   cuts.discr <- lapply(apply(rbind("nbcat"=categories, data), 2, function(x) { y <- x[1]; x <- x[-1]; return(list(quantile(x=x, probs=seq(0, 1, length.out=y+1), na.rm=TRUE)[-c(1, y+1)])) }), function(x) { return(x[[1]]) })
		   } else { cuts.discr <- categories }
## discretize the actual gene expression data using cutoffs identified from the training set
		   data <- data.discretize(data=data, cuts=cuts.discr)
		   nbcat <- sapply(cuts.discr, length) + 1
		   categories <- sapply(nbcat, function(x) { return(list(1:x))})
		   } else {
## data are already discretized
		   nbcat <- apply(data, 2, function(x) { return(length(sort(unique(x)))) })
		   categories <- lapply(apply(data, 2, function(x) { return(list(sort(unique(x)))) }), function(x) { return(x[[1]]) })
## not ideal solution since we may miss a category not present in the data
		   }
## data are discretized and nbcat is a list with the corresponding categories
		   
		   bnet <- .fit.catnet(data=data, categories=categories, perturbations=perturbations, priors=priors, priors.weight=priors.weight, maxparents=maxparents, maxparents.push=maxparents.push, seed=seed, bayesnet.maxcomplexity=bayesnet.maxcomplexity, bayesnet.maxiter=bayesnet.maxiter)
		   edgerel <- bnet$edge.relevance
## remove edge relevance from the network object to avoid redunddancy
		   bnet <- bnet[!is.element(names(bnet), "edge.relevance")]
#		if(retoptions=="all") {
			return(list("method"=method, "topology"=t(cnMatParents(bnet$model)), "topology.coeff"=NULL, "net"=bnet, "edge.relevance"=edgerel))
#						} else {
#		   return(list("method"=method, "topology"=.bayesnet2topo(net=bnet), "topology.coeff"=NULL, "edge.relevance"=edgerel))
#						}
		   }
		   }, 
		   
		   "regrnet"={
		   if(ensemble){
		   if(sum(perturbations)>0){
				  stop("perturbations cannot be considered in the ensemble approach")
		   }
## fit an ensemble regression model
		   if(missing(predn)){predn<-seq(1,ncol(data))}else if(length(intersect(predn,colnames(data)))>0){predn<-match(predn,colnames(data))}
		   rep_boot <- 200
		   if(ensemble.model=="best"){
		   vec_ensemble <- .Call("mrmr_ensemble", data.matrix(data),as.integer(is.na(data)),maxparents, ncol(data), nrow(data), predn, length(predn), rep_boot, ensemble.maxnsol, -1000)
		   }else if (ensemble.model=="full"){
		   vec_ensemble <- .Call("mrmr_ensemble_remove", data.matrix(data),as.integer(is.na(data)),maxparents, ncol(data), nrow(data), predn, length(predn), rep_boot, ensemble.maxnsol, -1000)
		   }
#net <- .extract.adjacency.ensemble(data,vec_ensemble,predn)
		   models.equiv <- .extract.all.parents(data,vec_ensemble,maxparents,predn)
		   if(causal){
		   res.causal<- .rank.genes.causal.ensemble(models.equiv,data)
		   }else{
		   res.causal<- .rank.genes.ensemble(models.equiv,data)
		   }
		   res.regrnet.ensemble<- .fit.regrnet.causal.ensemble(res.causal,models.equiv,data,priors=priors,priors.weight=priors.weight,maxparents=maxparents)
		   
		   tmp.res<-list("method"="regrnet", "topology"=.regrnet2topo.ensemble(net=res.regrnet.ensemble,coefficients=FALSE), "topology.coeff"=.regrnet2topo.ensemble(net=res.regrnet.ensemble,coefficients=TRUE), "edge.relevance"=res.regrnet.ensemble$edge.relevance)

			return(.list2matrixens(tmp.res))
		   } else {
## fit regression model
		   if(priors.count) {
## scale the priors
## truncate the largest counts
		   pp <- quantile(x=abs(priors[priors != 0]), probs=0.95, na.rm=TRUE)
		   if(length(pp) > 0) {
## truncate the distribution
		   priors[priors > pp] <- pp
		   priors[priors < -pp] <- -pp
		   }
## rescale the priors to have values in [-1, 1]
		   priors <- priors / max(abs(priors), na.rm=TRUE)
		   } else {
		   if(max(priors, na.rm=TRUE) > 1 || min(priors, na.rm=TRUE) < 0) { stop("'priors' should contain probabilities of interactions if 'priors.count' is FALSE!") }
## priors are probabilties that should be rescaled in [-1, 1]
		   priors <- (priors - 0.5) * 2
		   }
		   bnet <- .fit.regrnet.causal(data=data, perturbations=perturbations, priors=priors, predn=predn, maxparents=maxparents, maxparents.push=maxparents.push, priors.weight=priors.weight, causal=causal, regrmodel=regrmodel, seed=seed)
#return(bnet)
			## topology
			ttopo <- .regrnet2topo(net=bnet, coefficients=FALSE)
			## rescale score for each edge to [0, 1] 
			edgerel <- (bnet$edge.relevance + 1) / 2
			## absent edges are assigned a relevance score of 0
			edgerel[ttopo != 1] <- 0
			## remove edge relevance from the network object to avoid redundancy
		   bnet <- bnet[!is.element(names(bnet), "edge.relevance")]
#			if(retoptions=="all") {
#	return(list("method"=method, "topology"=.regrnet2topo(net=bnet, coefficients=FALSE), "topology.coeff"=.regrnet2topo(net=bnet, coefficients=TRUE), "net"=bnet, "edge.relevance"=edgerel))
#							} else {
			return(list("method"=method, "topology"=ttopo, "topology.coeff"=.regrnet2topo(net=bnet, coefficients=TRUE), "edge.relevance"=edgerel))
#							}
		   }
		   }
		   )
}
