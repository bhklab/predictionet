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
## seed: seed to make the function fully deterministic
## bayesnet.maxcomplexity: maximum complexity for bayesian network inference
## bayesnet.maxiter: maximum number of iterations for bayesian network inference

`netinf` <- function(data, categories, perturbations, priors, predn, priors.count=TRUE, priors.weight=0.5, maxparents=3, subset, method=c("regrnet", "bayesnet"), ensemble=FALSE, ensemble.model=c("full", "best"), ensemble.maxnsol=3, causal=TRUE, seed, bayesnet.maxcomplexity=0, bayesnet.maxiter=100, verbose=FALSE) {
## select more genes to find a number of cauqal variables close or more than maxparents
	if(ncol(data) < 2) { stop("Number of variables is too small to infer a network!") }
	if(!missing(predn) && !is.null(predn) && length(predn) < 2) { stop("length of parameter 'predn' should be >= 2!") }
	maxparents <- ifelse(maxparents >= ncol(data), ncol(data)-1, maxparents)
	maxparents <- ifelse(maxparents < 2, 2, maxparents)
	method <- match.arg(method)
	ensemble.model <- match.arg(ensemble.model)
	if(missing(perturbations) || is.null(perturbations)) {
		perturbations <- matrix(FALSE, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	} else {
		if(nrow(perturbations) == 1) {
			perturbations[1, ] <- as.logical(perturbations[1, ])
		} else { perturbations <- apply(perturbations, 2, as.logical) }
		dimnames(perturbations) <- dimnames(data)
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
		   if(ensemble) {
		   stop("ensemble bayesian network inference is not implemented yet!")
		   } else {
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
		   
			   bnet <- .fit.catnet(data=data, categories=categories, perturbations=perturbations, priors=priors, priors.weight=priors.weight, maxparents=maxparents, seed=seed, bayesnet.maxcomplexity=bayesnet.maxcomplexity, bayesnet.maxiter=bayesnet.maxiter)
			   edgerel <- bnet$edge.relevance
				return(list("method"=method, "ensemble"=ensemble, "topology"=t(cnMatParents(bnet$model)), "edge.relevance"=edgerel, "edge.relevance.global"=NULL))
			}
		  }, 
		   
		   "regrnet"={
		   if(ensemble){
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
		   
				if(sum(perturbations)>0){
					stop("perturbations cannot be considered in the ensemble approach")
				}
## fit an ensemble regression model
				if(missing(predn)){predn <- seq(1,ncol(data))}else if(length(intersect(predn,colnames(data)))>0){predn <- match(predn,colnames(data))}
				rep_boot <- 200
				if(ensemble.model=="best"){
					vec_ensemble <- .Call("mrmr_ensemble", data.matrix(data),as.integer(is.na(data)),maxparents, ncol(data), nrow(data), predn, length(predn), rep_boot, ensemble.maxnsol, -1000)
				}else if(ensemble.model=="full"){
		   vec_ensemble <- NULL
		   models.equiv <- NULL
		   
		   for(ind in 1:length(predn)){
		   if(verbose){
				print(paste("model for node",predn[ind],"is being built!"))
		   }
		   vec_ensemble <- c(vec_ensemble,list(.Call("mrmr_ensemble_remove", data.matrix(data),as.integer(is.na(data)),maxparents, ncol(data), nrow(data), predn[ind], 1, rep_boot, ensemble.maxnsol, -1000)))
		   models.equiv <- cbind(models.equiv,.extract.all.parents(data,vec_ensemble[[ind]],maxparents,predn[ind]))
		   }
		  
		   }
				if(verbose) { message("ensemble done") }
#	models.equiv <- .extract.all.parents(data,vec_ensemble,maxparents,predn)
				if(causal){
					res.causal <- .rank.genes.causal.ensemble(models.equiv,data)
				}else{
					res.causal <- .rank.genes.ensemble(models.equiv,data)
				}
				res.regrnet.ensemble <- .fit.regrnet.causal.ensemble(res.causal,models.equiv,data,priors=priors,priors.weight=priors.weight,maxparents=maxparents)
			## identify topology
			topores <- matrix(0, nrow=length(res.regrnet.ensemble$varnames), ncol=length(res.regrnet.ensemble$input), dimnames=list(res.regrnet.ensemble$varnames, res.regrnet.ensemble$input))
			for(i in 1:length(res.regrnet.ensemble$input)){
				topores[,i] <- res.regrnet.ensemble$edge.relevance[[i]]
			}
			topores[topores != 0] <- 1
			## build network object
		   
		   tmp.res <- list("method"="regrnet", "ensemble"=ensemble, "topology"=topores, "edge.relevance"=res.regrnet.ensemble$edge.relevance, "edge.relevance.global"=NULL)
		   tmp.res <- .list2matrixens(tmp.res)
		   
		   ind.rebuild <- NULL
		   unique.colnames <- unique(colnames(tmp.res$topology))
		   
		   for(i in 1:length(unique.colnames)){
				ind <- which(colnames(tmp.res$topology)==unique.colnames[i])
				if(length(ind)==1){
					ind.rebuild <- c(ind.rebuild,ind)
				}else{
					tmp.topo <- tmp.res$topology[,ind,drop=FALSE]
					ind.mult <- which(duplicated(tmp.topo, MARGIN = 2))
					if(length(ind.mult)>0){
						ind.rebuild <- c(ind.rebuild,ind[-ind.mult])
					}else{
						ind.rebuild <- c(ind.rebuild,ind)
					}
				}
		   }
		   tmp.res$topology <- tmp.res$topology[,ind.rebuild,drop=FALSE]
		   tmp.res$edge.relevance <- tmp.res$edge.relevance[,ind.rebuild,drop=FALSE]

		
		   return(tmp.res)

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
		   bnet <- .fit.regrnet.causal2(data=data, perturbations=perturbations, priors=priors, predn=predn, maxparents=maxparents, priors.weight=priors.weight, causal=causal, seed=seed)

## topology
			topores <- matrix(0, nrow=length(bnet$varnames), ncol=length(bnet$input), dimnames=list(bnet$varnames, names(bnet$input)))
			for(i in 1:length(bnet$input)) { topores[bnet$input[[i]], names(bnet$input)[i]] <- 1 }
			## rescale score for each edge to [0, 1] 
		   score <- bnet$edge.relevance
		   
			edgerel <- (bnet$edge.relevance + 1) / 2
			## absent edges are assigned a relevance score of 0
			edgerel[topores != 1] <- 0
		   return(list("method"=method, "ensemble"=ensemble, "topology"=topores, "edge.relevance"=edgerel, "edge.relevance.global"=score))
		   }
		})
}
