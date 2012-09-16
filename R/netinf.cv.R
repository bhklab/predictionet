### Function inferring a network
## data: matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## categories: Either a single integer or a vector of integers specifying the number of categories used to discretize each variable (data are then discretized using equal-frequency bins) or a list of cutoffs to use to discretize each of the variables in 'data' matrix. If method='bayesnet', this parameter should be specified by the user.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## priors: matrix of prior information avalable for gene-gene interaction (parents in columns, children in rows). Values may be probabilities or any other weights (citations count for instance)
## priors.count: 'TRUE' if priors specified by the user are number of citations for each interaction, 'FALSE' if probabilities are reported instead
## priors.weight: real value in [0, 1] specifying the weight to put on the priors (0=only the data are used, 1=only the priors are used to infer the topology of the network). If 'priors.weight' is missing it will be optimized gene by hene in an automatic way.
## maxparents: maximum number of parents allowed for each gene
## subset: vector of indices to select only subset of the observations
## method: "bayesnet" for bayesian network inference with the catnet package, "regrnet" for regression-based network inference
## nfold: number of folds for the cross-validation
## causal: 'TRUE' if the causality should be inferred from the data, 'FALSE' otherwise }
## seed: set the seed to make the cross-validation and network inference deterministic
`netinf.cv` <- 
function(data, categories, perturbations, priors, predn, priors.count=TRUE, priors.weight=0.5, maxparents=3, subset, method=c("regrnet", "bayesnet"),ensemble=FALSE, ensemble.maxnsol=3, predmodel=c("linear", "linear.penalized", "cpt"), nfold=10, causal=TRUE, seed, bayesnet.maxcomplexity=0, bayesnet.maxiter=100, verbose=FALSE) {
	if(!missing(seed)) { set.seed(seed) }
	method <- match.arg(method)
	predmodel <- match.arg(predmodel)
	if(missing(perturbations) || is.null(perturbations)) {
		perturbations <- matrix(FALSE, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	} else {
		if(nrow(perturbations) == 1) {
			perturbations[1, ] <- as.logical(perturbations[1, ])
		} else { perturbations <- apply(perturbations, 2, as.logical) }
		dimnames(perturbations) <- dimnames(data)
	}	
	if(!missing(subset)) {
		## select subset of the data (observations)
		data <- data[subset, , drop=FALSE]
		perturbations <- perturbations[subset, , drop=FALSE]
	}
	if(missing(predn) || is.null(predn) || length(predn) == 0) { predn <- 1:ncol(data) } else { if(is.character(predn)) { predn <- match(predn, dimnames(data)[[2]]) } else { if(!is.numeric(predn) || !all(predn %in% 1:ncol(data))) { stop("parameter 'predn' should contain either the names or the indices of the variables to predict!")} } }
	nn <- dimnames(data)[[2]][predn] ## variables (genes) to fit during network inference
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
		## discretize the actual gene expression data
		data.discr <- data.discretize(data=data, cuts=cuts.discr)
		categories <- cuts.discr
	} else {
		data.discr <- data
		categories <- lapply(apply(data, 2, function(x) { return(list(sort(unique(x)))) }), function(x) { return(x[[1]]) })
		names(categories) <- colnames(data)
	}
	
	## infer network from the whole dataset
	mynetglobal <- netinf2(data=data, categories=categories, perturbations=perturbations, priors=priors, predn=predn, priors.count=priors.count, priors.weight=priors.weight, maxparents=maxparents, method=method, ensemble=ensemble, causal=causal)
	## and the global topology
	mytopoglobal <- mynetglobal$topology

	## compute folds for cross-validation
	if (nfold == 1) {
		k <- 1
		nfold <- nrow(data)
	} else { k <- floor(nrow(data) / nfold) } 
	## randomize observations to prevent potential bias in the cross-validation
	smpl <- sample(nrow(data))
	
	## perform cross-validation
	mymcc <- mynrmse <- myr2 <- mytopo <- mynets <- NULL
	mytopo2 <- NULL
	edge.relevance.cv <- NULL
	
	for (i in 1:nfold) {
		if(ensemble & verbose) {
			message(paste("fold: ",i," out of ",nfold,"folds"))
		}
		## fold of cross-validation
		if (i == nfold) { s.ix <- smpl[c(((i - 1) * k + 1):nrow(data))] }	else { s.ix <- smpl[c(((i - 1) * k + 1):(i * k))] }
		## s.ix contains the indices of the test set
		
		## infer network from training data and priors

		mynet <- netinf2(data=data[-s.ix, , drop=FALSE], categories=categories, perturbations=perturbations[-s.ix, , drop=FALSE], priors=priors, predn=predn, priors.count=priors.count, priors.weight=priors.weight, maxparents=maxparents, method=method, ensemble=ensemble, ensemble.maxnsol=ensemble.maxnsol, causal=causal, bayesnet.maxcomplexity=bayesnet.maxcomplexity, bayesnet.maxiter=bayesnet.maxiter)
		mynet <- net2pred(net=mynet, data=data[-s.ix, , drop=FALSE], categories=categories, predn=predn, perturbations=perturbations[-s.ix, , drop=FALSE], method=predmodel)
		mynets <- c(mynets, list(mynet))
		
		#mynet <- c(mynet,list("topology.coeff"=.regrnet2matrixtopo(net=mynet)))
		mytopo2 <- c(mytopo2, list(mynet$topology))

		## edge relevance score
		edge.relevance.cv <- c(edge.relevance.cv,list(mynet$edge.relevance))
		## compute predictions
		mynet.pred <- netinf.predict(net=mynet, data=data[s.ix, , drop=FALSE], categories=categories, perturbations=perturbations[s.ix, , drop=FALSE], predn=predn, method=predmodel)
		## combine the multiple predictions computed by the ensemble method by averaging
		uu <- !duplicated(colnames(mynet.pred))
		if(sum(uu) < ncol(mynet.pred)) {
			ttt <- matrix(NA, nrow=nrow(data[s.ix, , drop=FALSE]), ncol=ncol(data[s.ix, , drop=FALSE]), dimnames=dimnames(data[s.ix, , drop=FALSE]))
			ttt[ , colnames(mynet.pred)[uu]] <- mynet.pred[ , uu]
			uu <- sort(unique(colnames(mynet.pred)[duplicated(colnames(mynet.pred))]))
			for(u in 1:length(uu)) {
				ttt[ , uu[u]] <- apply(mynet.pred[ , is.element(colnames(mynet.pred), uu[u])], 1, mean, na.rm=TRUE)
			}
			mynet.pred <- ttt
		}
		resnull <- rep(NA, ncol(data))
		names(resnull) <- colnames(data)
		
		## performance estimation: R2
		mynet.r2 <- resnull
		if(predmodel %in% c("linear", "linear.penalized")) {
			mynet.r2 <- pred.score(data=data[s.ix, , drop=FALSE], pred=mynet.pred, method="r2")
		}
			
		## performance estimation: NRMSE
		mynet.nrmse <- resnull
		if(predmodel %in% c("linear", "linear.penalized")) {
			mynet.nrmse <- pred.score(data=data[s.ix, , drop=FALSE], pred=mynet.pred, method="nrmse")
		}
		## performance estimation: MCC
		pred.discr <- mynet.pred
		if(predmodel %in% c("linear", "linear.penalized") && !missing(categories)) {
			## discretize predicted gene expression data 
			pred.discr <- data.discretize(data=mynet.pred, cuts=cuts.discr[dimnames(mynet.pred)[[2]]])
		}
		if(predmodel %in% c("cpt") || (predmodel %in% c("linear", "linear.penalized") && !missing(categories))) { 
			mynet.mcc <- pred.score(data=data.discr[s.ix, dimnames(mynet.pred)[[2]], drop=FALSE], pred=pred.discr, method="mcc")
		} else { 
			mynet.mcc <- rep(NA, length(nn))
			names(mynet.mcc) <- nn
		}
		
		## adjacency matrix
		topol <- mynet$topology
		## save results
		myr2 <- rbind(myr2, mynet.r2)
		mynrmse <- rbind(mynrmse, mynet.nrmse)
		mymcc <- rbind(mymcc, mynet.mcc)
		mytopo <- c(mytopo, list(topol))
	}
	dimnames(myr2)[[1]] <- dimnames(mynrmse)[[1]] <- dimnames(mymcc)[[1]] <- names(mytopo) <- names(mynets) <- names(edge.relevance.cv) <- paste("fold", 1:nfold, sep=".")
	
	## compute stability for each edge
	edgestab <- edgestab2 <- matrix(0, nrow=nrow(mytopo[[1]]), ncol=ncol(mytopo[[1]]), dimnames=dimnames(mytopo[[1]]))

	if(ensemble){
			edgestab <- .topo2stab(mytopo)
			edgestab2 <- .stab.cv2stab(edgestab,mytopoglobal)
	}else{
		for(i in 1:length(mytopo)) {
			edgestab <- edgestab + mytopo[[i]]
		}
		edgestab <- edgestab / length(mytopo)
		edgestab2[mytopoglobal == 1] <- edgestab[mytopoglobal == 1]

	}
	## report stability of edges present in the global network
	names(mytopo2) <- paste("fold", 1:nfold, sep=".")
	return(list("method"=method, "topology"=mytopoglobal, "topology.cv"=mytopo, "prediction.score.cv"=list("r2"=(myr2), "nrmse"=(mynrmse), "mcc"=(mymcc)), "edge.stability"=edgestab2, "edge.stability.cv"=edgestab, "edge.relevance"=mynetglobal$edge.relevance, "edge.relevance.cv"=edge.relevance.cv))
}
