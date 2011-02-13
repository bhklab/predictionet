## unctions of our netinf package

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
`netinf` <- 
function(data, categories, perturbations, priors, predn, priors.count=TRUE, priors.weight=0.5, maxparents=3, subset, method=c("regrnet", "bayesnet"), causal=TRUE, seed) {
	
	if(!missing(predn) && !is.null(predn) && length(predn) < 2) { stop("length of parameter 'predn' should be >= 2!") }
	if(causal && maxparents > (ncol(data) * 0.5)) { warning("Maximum number of parents may be too large, causal inference requires sparsity in the inferred network; please decrease maxparents parameter for better results!") }
	method <- match.arg(method)
	if(missing(perturbations)) {
		## create matrix of no perturbations
		perturbations <- matrix(0, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	}
	if(missing(priors)) {
		## create a matrix of priors count = 0
		priors <- matrix(0, nrow=ncol(data), ncol=ncol(data), dimnames=list(dimnames(data)[[2]], dimnames(data)[[2]]))
	}
	if(priors.weight < 0 || priors.weight > 1) { stop("'priors.weight' should be in [0, 1]!") }
	if(!missing(subset)) {
		data <- data[subset, , drop=FALSE]
		perturbations <- perturbations[subset, , drop=FALSE]
	}
	switch(method, 
	"bayesnet"={
		cat("Bayesian network inference is not implemented yet!\n")
		return(0)
	}, 
	"regrnet"={
		## fit regression model
		if(priors.count) {
			## scale the priors
			## truncate the largest counts
			pp <- quantile(x=abs(priors[priors != 0]), probs=0.99, na.rm=TRUE)
			if(length(pp) > 0) {
				priors[priors > pp] <- pp
				priors[priors < -pp] <- -pp
			}
			## rescale the priors to have values in [-1, 1]
			priors <- priors / max(abs(priors), na.rm=TRUE)
		} else if(max(priors)>1 || min(priors)< -1) {
			stop("'priors' should contain count values (e.g. number of citations) for method 'regrnet' if 'priors.count' is TRUE and values in [-1,1] if'priors.count' is FALSE!")
		}
		bnet.regr <- fit.regrnet.causal(data=data, perturbations=perturbations, priors=priors, predn=predn, maxparents=maxparents, priors.weight=priors.weight, causal=causal, seed=seed)
		return(list("method"=method, "net"=bnet.regr))
	})
			
}

### Function predicting values of nodes from a previously inferred network
`netinf.predict` <- 
function(net, data, categories, perturbations, subset, predn) {
	switch(net$method, 
		"bayesnet"={
			cat("Bayesian network inference is not implemented yet!\n")
			return(0)
		}, 
		"regrnet"={
		if(missing(categories)){
			res <- pred.onegene.regrnet.fs(net=net$net, data=data, perturbations=perturbations, subset=subset, predn=predn)
			return(res)
		} else {
			return(0)
		}
	})
}

### function fitting a regression model for each gene in the data
## data: Matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## priors: matrix of prior information avalable for gene-gene interaction (parents in columns, children in rows). Values may be probabilities or any other weights (citations count for instance)
## predn: indices or names of genes to fit during network inference. If missing, all the genes will be used for network inference
## priors.weight: real value in [0, 1] specifying the weight to put on the priors (0=only the data are used, 1=only the priors are used to infer the topology of the network). If 'priors.weight' is missing it will be optimized gene by hene in an automatic way.
## maxparents: maximum number of parents allowed for each gene
### returns a regression network 
`fit.regrnet.causal` <- 
function(data, perturbations, priors, predn, maxparents, priors.weight, causal=TRUE, seed) {
	require(penalized)
	
	if(!missing(seed)) { set.seed(seed) }
	
	if(missing(predn) || is.null(predn)) { predn <- 1:ncol(data) } else { if(is.character(predn)) { predn <- match(predn, dimnames(data)[[2]]) } else { if(!is.numeric(predn) || !all(predn %in% 1:ncol(data))) { stop("parameter 'predn' should contain either the names or the indices of the variables to predict!")} } }
	nn <- dimnames(data)[[2]][predn] ## variables (genes) to fit during network inference
	
	########################
	### data preparation
	########################
	dd <- data.frame(data)
	diag(priors) <- 0 ## no self-loops
	
	regrnet <- NULL
	
	########################
	### find ranked parents for all genes
	########################
	mat.ranked.parents <- rank.genes.causal.perturbations(priors=priors, data=dd, perturbations=perturbations, predn=predn, priors.weight=priors.weight, maxparents=maxparents, causal=causal)
	myparents <- lapply(apply(mat.ranked.parents, 1, function(x) { x <- x[!is.na(x)]; if(length(x) == 0) { x <- NULL };  return(list(x)); }), function(x) { return(x[[1]]) })

	########################
	### build regression model for each gene
	########################
	for(i in 1:length(nn)) {
		### consider only those genes for the model which are ranked (not the NA entries)
		modelv <- myparents[[nn[i]]]
		if(is.null(modelv)) { ## none variables have been selected
			mm.reg <- mean(dd[ , nn[i]], na.rm=TRUE) ## null model (only intercept)
		} else {
			modelv <- sort(modelv) ## sort variables by name and remove NA
			#ff <- formula(sprintf("%s ~ %s%s", nn[i], ifelse(is.null(modelv), "1", "1 + "), paste(modelv, collapse=" + ")))
			ff <- formula(sprintf("%s ~ 1 + %s", nn[i], paste(modelv, collapse=" + ")))
			optlambda1 <- optL1(response=ff, data=dd, model="linear", lambda2=0, minlambda1=1, maxlambda1=10, trace=FALSE, fold=10)
			mm.reg <- penalized(response=ff, data=dd, model="linear", lambda1=optlambda1$lambda, lambda2=0, trace=FALSE)
			## cor(data[ , nn[i]], predict(object=mm.reg, data=dd)[ , 1])
			## compare with lm
			## mm <- lm(formula=ff, data=dd)
			## cor(data[ , nn[i]], predict(mm, newdata=dd))
		}
		regrnet <- c(regrnet, list(mm.reg))
	}
	names(regrnet) <- nn
	
	return(list("varnames"=colnames(data), "input"=myparents, "model"=regrnet))
}

### function ranking the genes based on mrmr+prior+causality for each gene as a target
## data: Matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## priors: matrix of prior information avalable for gene-gene interaction (parents in columns, children in rows). Values may be probabilities or any other weights (citations count for instance)
## predn: indices of genes to fit during network inference. If missing, all the genes will be used for network inference
## priors.weight: real value in [0, 1] specifying the weight to put on the priors (0=only the data are used, 1=only the priors are used to infer the topology of the network). If 'priors.weight' is missing it will be optimized gene by hene in an automatic way.
## maxparents: maximum number of parents allowed for each gene
### returns matrix: each row corresponds to a target gene, columns contain the ranked genes (non NA values; gene names)
`rank.genes.causal.perturbations` <- 
function(priors, data, perturbations, predn, priors.weight, maxparents, causal=TRUE) {

	########################
	### initialize variables
	########################
	res.netw <- matrix(0, nrow=ncol(data), ncol=length(predn), dimnames=list(dimnames(data)[[2]], dimnames(data)[[2]][predn]))
	res <- matrix(NA, nrow=length(predn), ncol=ncol(data) - 1, dimnames=list(dimnames(data)[[2]][predn], seq(1, ncol(data) - 1)))

	data.original <- data
	
	if(all(apply(perturbations, 2, sum, na.rm=TRUE) == 0) & sum(is.na(perturbations)) == 0) {
		## no perturbations
		## estimate the mutual information matrix for all pairs of genes
		mim.global <- build2.mim(data.original)
		## compute mrmr score for the target genes
		mrmr.global <- .Call("mrnet_adapted", mim.global, ncol(data.original), -1000)
		mrmr.global <- matrix(mrmr.global, nrow=nrow(mim.global), ncol=ncol(mim.global))
		dimnames(mrmr.global) <- list(colnames(data.original), colnames(data.original))
	}
	
	for(i in 1:length(predn)) {
		data <- data.frame(data.original)
		
		## remove observation where the target variabl ewas perturbed (no link between observation and target variable)
		ind <- which(is.na(perturbations[, predn[i]]) | perturbations[, predn[i]] == 1)
		if(length(ind)>0) {
			data <- data[-ind,]
			## since the target gene has been perturbed in some of the experiments we should recompute the mim and mrmr matrices
			## estimate the mutual information matrix for all pairs of genes
			mim <- build2.mim(data)
			## compute mrmr score for the target genes
			mrmr <- .Call("mrnet_adapted", mim, ncol(data), -1000)
			mrmr <- matrix(mrmr, nrow=nrow(mim), ncol=ncol(mim))
			dimnames(mrmr) <- list(colnames(data),colnames(data))
		} else {
			mim <- mim.global
			mrmr <- mrmr.global
		}

		########################
		### compute mrmr score matrix
		########################
		
		
		## normalization of mrmr scores
		mrmr <- mrmr/max(mrmr)
		## make score symmetric
		for(k in 1:ncol(data)) {
			for(j in k:ncol(data)) {
				if(mrmr[k,j]>0 && mrmr[k,j]>mrmr[j,k]) {
					mrmr[j,k] <- mrmr[k,j]
				} else if(mrmr[j,k]>0 && mrmr[k,j]<mrmr[j,k]) {
					mrmr[k,j] <- mrmr[j,k]
				} else if(mrmr[k,j]<0 & mrmr[j,k]<0) {
					if(mrmr[k,j]>mrmr[j,k]) {
						mrmr[j,k] <- mrmr[k,j]
					} else {
						mrmr[k,j] <- mrmr[j,k]
					}
				}
			}
		}
		## avoid mrmr score of exactly 0
		mrmr[mrmr == 0] <- 1e-16
		## mrmr scores on the diagonal should be equal to 0
		diag(mrmr) <- 0
		score <- matrix(0, nrow=ncol(data), ncol=ncol(data), dimnames=list(dimnames(data)[[2]], dimnames(data)[[2]]))
		### limit the number of adjacent variables to maxparents
		for(j in 1:ncol(data)){
			tmp.s <- order(mrmr[ ,j],decreasing=TRUE)
			## we do not want to select the target gene
			tmp.s <- tmp.s[-which(tmp.s == j)]
			score[tmp.s[1:maxparents],j] <- mrmr[tmp.s[1:maxparents],j]
		}

		if(causal) {
			########################
			## find all adjacencies of node to predict in score
			########################
			tmp.adj <- matrix(0,nc=ncol(score),nr=nrow(score))

			tmp.adj[predn[i],] <- tmp.adj[,predn[i]] <- score[,predn[i]]

			########################
			## determine all triplets X.k-X.j-X.l from these adjecencies
			########################
			trip <- (network2triplets(tmp.adj))
			if(length(trip)>0){
				if(length(trip)==3){
					trip <- (t(as.matrix(trip)))
				}
				########################
				## compute I(X.k,X.l)-I(X.k,X.l|X.j)
				########################
		
				trip <- as.data.frame(cbind(trip,(get.ii4triplets.gaussian(data,trip))))
				########################
				## remove those triplets for which the outer nodes are connected
				########################
				ind.rm <- NULL
				for(j in 1:nrow(trip)){
					if(score[trip[j,1],trip[j,3]] !=0 || score[trip[j,3],trip[j,1]] !=0) {
						ind.rm <- c(ind.rm,j)
					}
				}
				if(!is.null(ind.rm)) { trip <- (trip[-ind.rm, , drop=FALSE]) }
				########################
				## if there are triplets for which I(X,Z)-I(X,Z|Y)<0, add these to the list of parents
				########################
				parents <- NULL
				if(nrow(trip)>0){
					for(j in 1:nrow(trip)){
						if(trip[j,4]<0){
							parents <- c(parents,trip[j,1],trip[j,3])
						}
					}
					parents <- unique(parents) 
				}
			}
		} else {
			## no causality inference, all connected are parents
			parents <- which(score[ , predn[i]] != 0)
			if(length(parents) == 0) { parents <- NULL }
		}
		# parents into network matrix
		res.netw[parents,i] <- score[parents,predn[i]]
	}
	
	########################
	## w*M+(1-w)*P weighted score from prior knowledge and data
	########################
	score <- (1-priors.weight)*res.netw+priors.weight*priors[ , predn, drop=FALSE]
	## no self-loops
	diag(score[dimnames(score)[[2]], dimnames(score)[[2]]]) <- 0
	for(j in 1:length(predn)){
		tmp.s <- sort(score[,j],decreasing=TRUE,index.return=TRUE)
		ind.rm <- (which(tmp.s$x<=0))
		if(length(tmp.s$x[-ind.rm])>0){
			res[j,(1:(length(tmp.s$x[-ind.rm])))] <- names(tmp.s$x[-ind.rm])
			if(length(tmp.s$x[-ind.rm])>maxparents){
				res[j,((maxparents+1):(length(tmp.s$x[-ind.rm])))] <- NA
			}
		}
	}
	return(res)
}

### function predicting one gene at the time for the test data using the regression model built on the training data
## net: regression model for each gene
## data: Matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## predn: indices or names of variables (genes) to predict
### return continuous predictions for each gene
`pred.onegene.regrnet.fs` <- 
function(net, data, perturbations, subset, predn) {
	if(missing(perturbations) || is.null(perturbations)) {
		perturbations <- matrix(FALSE, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	} else {
		perturbations <- apply(perturbations, 2, as.logical)
		dimnames(perturbations) <- dimnames(data)
	}
	if(!missing(subset)) {
		data <- data[subset, , drop=FALSE]
		perturbations <- perturbations[subset, , drop=FALSE]
	}
	if(missing(predn) || is.null(predn)) { predn <- match(names(net$input), dimnames(data)[[2]]) } else { if(is.character(predn)) { predn <- match(predn, dimnames(data)[[2]]) } else { if(!is.numeric(predn) || !all(predn %in% 1:ncol(data))) { stop("parameter 'predn' should contain either the names or the indices of the variables to predict!")} } }
	if(!all(is.element(predn, match(names(net$input), dimnames(data)[[2]])))) { stop("some genes cannot be predicted because they have not been fitted in the network!")}
	nnix <- dimnames(data)[[2]][predn] ## variables to predict
	
	## matrix to store the predictions
	preds <- matrix(NA, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	for (i in 1:length(nnix)) {
		## extract the fitted model
		model.i <- net$model[[nnix[i]]]
		if(!is.numeric(model.i)) {
			## non null model
			preds[ , nnix[i]] <- predict(object=model.i, data=data.frame(data))[ , 1]
		} else { preds[ , nnix[i]] <- model.i }
	}
	preds[perturbations] <- NA
	return(preds)
}

### function transforming a regression model into an adjacency matrix
## net: regression model for each gene
## coeffcients: if TRUE, the coefficients are extracted from the regression models for each existing edges,  boolean representing the presence of the edges otherwise
### returns adjacency matrix: parents in rows and children in columns (res[i,j]==1 means edge from i to j)
`regrnet2topo` <- 
function(net, coefficients=FALSE) {
	geneid <- names(net$input)
	nr <- length(geneid)
	net.model <- net$model
	res <- matrix(0, nrow=length(net$varnames), ncol=nr, dimnames=list(net$varnames,geneid))
	for (i in 1:nr) {
		model.i <- net.model[[i]]
		if(!is.numeric(model.i)) { ## non null model
			if(coefficients) {
				beta <- coefficients(object=model.i)
				res[names(beta)[-1], geneid[i]] <- beta[-1]
			} else { res[net$input[[i]], geneid[i]] <- 1 }
		}
	}
	return(res)
}

### function transforming a regression model into an adjacency matrix
## net: regression model for each gene
## coeffcients: if TRUE, the coefficients are extracted from the regression models for each existing edges,  boolean representing the presence of the edges otherwise
### returns adjacency matrix: parents in rows and children in columns (res[i,j]==1 means edge from i to j)
`bayesnet2topo` <- 
function(net) {
	stop("not implemented yet!")
}

### function transforming a regression model into an adjacency matrix
## net: regression model for each gene
## coeffcients: if TRUE, the coefficients are extracted from the regression models for each existing edges,  boolean representing the presence of the edges otherwise
### returns adjacency matrix: parents in rows and children in columns (res[i,j]==1 means edge from i to j)
`net2topo` <- 
function(net, ...) {
	switch(net$method, 
	"regrnet"={
		res <- regrnet2topo(net=net$net, ...)
	},
	"bayesnet"={
		res <- bayesnet2topo(net=net$net, ...)
	}, 
	{
		stop("method used to infer the network is unknown!")
	})
	return(res)
}

### function computing perfomance of prediction; methods include rnmse and mcc
### user has to provide data/pred accordingly (either both continuous or both discretized)
## data: matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## pred: Matrix of continuous or categorical values: predictions in rows, features in columns
## categories: list of categories for each of the nodes (genes) in 'data' matrix. Categories should be integers from 1 to n.
## method: scoring metric: nrmse, mcc,...
`pred.score` <- 
function(data, pred, categories, method=c("r2", "nrmse", "mcc")) {
	method <- match.arg(method)

	if(method %in% c("mcc")) {
		## need to discretize the data to compute this performance criterion
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
			pred <- data.discretize(data=pred, cuts=cuts.discr)
			nbcat <- sapply(cuts.discr, length) + 1
		} else {
			nbcat <- apply(rbind(data, pred), 2, max)
			## not ideal solution since we may miss a higher class not present in the data
		}
	} else {
		nbcat <- rep(0, ncol(data))
		names(nbcat) <- colnames(data)
	}
	
	myfoo <- function(data, pred, nbcat, method) {
		## data and pred should be vectors
		if(all(!complete.cases(data, pred))) { return(NA) }
		switch(method,
		"r2"={
			mss <- sum((pred - mean(pred))^2)
			resids <- (data - pred)
			rss <- sum(resids^2)
			myperf <- mss/(mss + rss)
		}, 
		"nrmse"={
			myperf <- sqrt(mean((data-pred)^2, na.rm=TRUE)) / (max(pred, na.rm=TRUE) - min(pred, na.rm=TRUE))
		},
		"mcc"={ ## Matthews Correlation Coefficient, equivalent to the phi coeffcient (see jurman2010unifying)
			tt2 <- as.table(matrix(0, nrow=nbcat, ncol=nbcat, dimnames=list(1:nbcat, 1:nbcat)))
			tt <- table("PREDS"=pred, "TRUTH"=data)
			tt2[dimnames(tt)[[1]], dimnames(tt)[[2]]] <- tt
			myperf <- (mcc(ct=tt2, nbcat=nbcat) + 1) / 2
			#if(is.na(myperf)) { myperf <- 0 }
		})
		return(myperf)
	}
	
	if(is.vector(data) && is.vector(pred)) {
		if(length(data) != length(pred)) { stop("length of 'data' and 'pred' should be equal!") }
		myperf <- myfoo(data=data, pred=pred, nbcat=nbcat, method=method)
	} else {
		if(!all(dim(data) == dim(pred))) { stop("dimensions of 'data' and 'pred' should be equal!") }
		myperf <- apply(X=rbind(nbcat, data, pred), MARGIN=2, FUN=function(x, foo) { nbcat <- x[1]; ll <- (length(x)-1)/2; data <- x[2:(ll+1)]; pred <- x[(ll+2):((2*ll)+1)]; return(foo(data=data, pred=pred, nbcat=nbcat, method=method)); }, foo=myfoo)
	}
	
	return(myperf)
}

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

## function to discretize data based on user specified cutoffs
## warning: the function requires an equal number of categories for each variable
## data: matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## cuts: list of cutoffs for each variable
`data.discretize` <- 
function(data, cuts) {
	if(!is.list(cuts)) { stop("'cuts' should be a list!") }
	if(length(cuts) != ncol(data)) { stop("the length of 'cuts' should be equal to the number of variables (collumns in data)!") }
	if(!is.null(names(cuts))) { cuts <- cuts[dimnames(data)[[2]]] }
	## transform the list of cutoffs in a matrix
	maxx <- max(unlist(lapply(cuts, function(x) { return(length(x[!is.na(x)])) })), na.rm=TRUE)
	cuts2 <- matrix(NA, nrow=ncol(data), ncol=maxx, dimnames=list(dimnames(data)[[2]], paste("cutoff", 1:maxx, sep=".")))
	for(i in 1:length(cuts)) { if(all(is.na(cuts[[i]]))) { cuts2[i, ] <- NA } else { cuts2[i, 1:sum(!is.na(cuts[[i]]))] <- cuts[[i]][!is.na(cuts[[i]])] } }
	
	myf <- function(x, y) {
		## data
		xx <- x[1:(length(x)-y)]
		## cutoffs
		cc <- x[(length(x)-y+1):length(x)]
		cc <- sort(cc[!is.na(cc)])
		if(length(cc) == 0) { return(xx) }
		xx2 <- rep(1, length(xx))
		for(i in 1:length(cc)) { xx2[!is.na(xx) & xx > cc[i]] <- xx2[!is.na(xx) & xx > cc[i]] + 1 }
		return(xx2)
	}
	
	res <- apply(rbind(data, t(cuts2)), 2, myf, y=maxx)
	dimnames(res) <- dimnames(data)
	return(res)
}

## function to compute the Matthews Correlation Coefficient (MCC) in a classification framework
## ct: confusion table (matrix or table)
## nbcat: the actual number of categories for the truth and the predictions
## Example:
## rr <- sapply(1:1000, function(x, cc) { set.seed(x); ctt <- table(sample(1:cc, 100, replace=TRUE), sample(1:cc, 100, replace=TRUE)); return(mcc(ct=ctt, nbcat=cc)); }, cc=3)
## hist(rr)
`mcc` <- 
function(ct, nbcat=nrow(ct)) {
	if(nrow(ct) != ncol(ct)) { stop("the confusion table should be squared!") }
	if(!(sum(ct)==sum(diag(ct))) &&  (length(which(apply(ct, 1, sum) == 0)) == (nbcat-1) & ((length(which(apply(ct, 2, sum) == 0)) != (nbcat-1)) | (length(which(apply(ct, 2, sum) == 0)) == (nbcat-1)))) || (length(which(apply(ct, 2, sum) == 0)) == (nbcat-1) & ((length(which(apply(ct, 1, sum) == 0)) != (nbcat-1)) | (length(which(apply(ct, 1, sum) == 0)) == (nbcat-1)) & sum(diag(ct)) == 0))) { ct <- ct + matrix(1, nc=nbcat, nr=nbcat) } ### add element to categories if nbcat-1 predictive categories do not contain elements. Not in case where all are correct!
	myx <- matrix(TRUE, nrow=nrow(ct), ncol=ncol(ct))
	diag(myx) <- FALSE
	if(sum(ct[myx]) == 0) { return(1) }
	myperf <- 0
	for(k in 1:nbcat) {
		for(m in 1:nbcat) {
			for(l in 1:nbcat) {
				myperf <- myperf + ((ct[k, k] * ct[m, l]) - (ct[l, k] * ct[k, m]))
			}
		}
	}
	aa <- 0
	for(k in 1:nbcat) {
		cc <- 0
		for(l in 1:nbcat) { cc <- cc + ct[l, k] }
		dd <- 0
		for(f in 1:nbcat) {
			for(g in 1:nbcat) { if(f != k) { dd <- dd + ct[g, f] } }
		}
		aa <- aa + (cc * dd)
	}
	bb <- 0
	for(k in 1:nbcat) {
		cc <- 0
		for(l in 1:nbcat) { cc <- cc + ct[k, l] }
		dd <- 0
		for(f in 1:nbcat) {
			for(g in 1:nbcat) { if(f != k) { dd <- dd + ct[f, g] } }
		}
		bb <- bb + (cc * dd)
	}
	
	myperf <- myperf / (sqrt(aa) * sqrt(bb))
	return(myperf)
}

### function finding all triplets in symmetric matrix X of the form A-B-C and A not connected with C
### X adjacency matrix (symmetric)
`network2triplets` <- 
function(X) {
### step 1: find all pairs (matrix in which the first element has a smaller number than the second)
### step 2: identify the triplets using matrix of pairs
	if(sum(X)!=0){
		n <- dim(X)[1]
		m <- length(X[X!=0])
		
		step1 <- get.pairs(X,n,m)
			
		pairs <- step1$pairs		 	# matrix of dimension m times 2, ordered pairs
		index <- step1$index		 	# index vector containing on position i the row number
		
		triplets <- matrix(0,nrow=90*m,ncol=3)
		cnt <- 1
		
		for(j in 1:m){
			
			p <- pairs[j,2]
			p2 <- pairs[j,1]
			
			p1 <- index[p]
			q1 <- index[p+1]-1
			
			indq1 <- 1
			while(q1==(-1)){
				q1 <- index[p+indq1]-1
				indq1 <- indq1+1
			}
			for(k in (p1:q1)){
				q2 <- pairs[k,2]
				
				p3 <- p2
				while(index[p3+1]==0){
					p3 <- p3+1
				}
				if(q2 > p2 && !((q2)%in% pairs[(index[p2]):(index[p3+1]-1),2])){ 
					triplets[cnt,] <- c(pairs[j,],q2)
					cnt <- cnt+1
				}
			}
		}
		return(triplets[triplets[,1]!=0,])
	}else{
		return(NULL)
	}		
}

### function returning an ordered matrix containing all pairs A-B in the symmetrix matrix X, 
### returns vector, element i contains position of start node i-... in matrix of pairs
### input, symmetrix matrix X containing 0 and 1, number of rows/columns, number of pairs 
### (=number of elements equal to 1)
`get.pairs` <- 
function(X,n,sum) {


	index.start.node <- rep(0,n)
	pairs <- matrix(0,nc=2,nr=sum)
	triplets <- matrix(0,nrow=30*sum,ncol=3)
	
	cnt <- 1
	ind <- 0
	
	index.start.node[1] <- 0
	for(i in seq(1,n)){
		for(j in seq(1,n)){
			if(X[i,j]!=0){
				if(i>ind){
					ind <- i
					index.start.node[i] <- cnt
				}
				pairs[cnt,] <- c(i,j)
				cnt <- cnt+1
			}
		}
	}
	return(list(pairs=pairs,index=c(index.start.node,(sum+1))))
}

### function to compute the interaction information of triplets of variables
### returns vector containing estimated interaction information values
## data: Matrix of continuous values (gene expressions for example); observations in rows, features in columns.
## triplets: matrix of triplets for which interaction information should be computed: rows containing triplets
`get.ii4triplets.gaussian` <- 
function(data, triplets, estimator=c("pearson", "spearman", "kendall")) {
### par: 	parameters to return, if 0 then both vector of values and vector of indices
### 		if !=0 then only value vector
	estimator <- match.arg(estimator)
	n <- dim(data)[2]
	nt <- dim(triplets)[1]

	allcor <- cor(data,method=estimator)
	
	vec.ci.val <- rep(0,nt)	
	for (l in 1:(nt)){
		tmp <- sort(triplets[l,]) ## use sorted triplets for performance reasons on the estimator
		i <- tmp[1]
		j <- tmp[2]
		k <- tmp[3]

		vec.ci.val[l] <-  -1/2 *log10(((1-allcor[i,j]^2)*(1-allcor[i,k]^2)*(1-allcor[j,k]^2))/(1+2*allcor[i,j]*allcor[i,k]*allcor[j,k]-allcor[i,j]^2-allcor[i,k]^2-allcor[j,k]^2))
	}
	return(vec.ci.val)
}

### function to compute the matrix of mutual information (inspired from minet)
### returns a matrix containing the mutual information between all pairs of variables
## datase: Matrix of continuous values (gene expressions for example); observations in rows, features in columns.
## estimator: type of correlation coefficient used to compute the mutual information
`build2.mim` <- 
function(dataset, estimator=c("pearson", "spearman", "kendall")) {
	estimator <- match.arg(estimator)
	if( estimator=="pearson" || estimator=="spearman" || estimator=="kendall") {
		  mim <- cor(dataset,method=estimator,use="complete.obs")^2
		  diag(mim) <- 0
		  maxi <- 0.999999
		  mim[which(mim>maxi)] <- maxi
		  mim  <- -0.5*log(1-mim)
	} else { stop("unknown estimator") }
		  
	mim[mim<0] <- 0
	return(mim)
}

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
function(data, categories, perturbations, priors, predn, priors.count=TRUE, priors.weight=0.5, maxparents=3, subset, method=c("regrnet", "bayesnet"), nfold=10, causal=TRUE, seed) {
	if(!missing(seed)) { set.seed(seed) }
	if(missing(perturbations)) {
		## create matrix of no perturbations
		perturbations <- matrix(0, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	}	
	if(!missing(subset)) {
		## select subset of the data (observations)
		data <- data[subset, , drop=FALSE]
		perturbations <- perturbations[subset, , drop=FALSE]
	}
	if(missing(predn) || is.null(predn)) { predn <- 1:ncol(data) } else { if(is.character(predn)) { predn <- match(predn, dimnames(data)[[2]]) } else { if(!is.numeric(predn) || !all(predn %in% 1:ncol(data))) { stop("parameter 'predn' should contain either the names or the indices of the variables to predict!")} } }
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
		mydata.discr <- data.discretize(data=data, cuts=cuts.discr)
		categories <- cuts.discr
	} else { mydata.discr <- data }
	
	## infer network from the whole dataset
	mynetglobal <- netinf(data=data, categories=categories, perturbations=perturbations, priors=priors, predn=predn, priors.count=priors.count, priors.weight=priors.weight, maxparents=maxparents, method=method, causal=causal)
	## and the global topology
	mytopoglobal <- regrnet2topo(net=mynetglobal$net, coefficients=FALSE)
	
	## compute folds for cross-validation
	if (nfold == 1) {
		k <- 1
		nfold <- nrow(data)
	} else { k <- floor(nrow(data) / nfold) }
	## randomize observations to prevent potential bias in the cross-validation
	smpl <- sample(nrow(data))

	## perform cross-validation
	mymcc <- mynrmse <- myr2 <- mytopo <- mynets <- NULL
	for (i in 1:nfold) {
		## fold of cross-validation
		if (i == nfold) { s.ix <- smpl[c(((i - 1) * k + 1):nrow(data))] }	else { s.ix <- smpl[c(((i - 1) * k + 1):(i * k))] }
		## s.ix contains the indices of the test set

		## infer network from training data and priors
		mynet <- netinf(data=data[-s.ix, , drop=FALSE], categories=categories, perturbations=perturbations[-s.ix, , drop=FALSE], priors=priors, predn=predn, priors.count=priors.count, priors.weight=priors.weight, maxparents=maxparents, method=method, causal=causal)
		mynets <- c(mynets, list(mynet))
		
		## compute predictions
		mynet.pred <-  netinf.predict(net=mynet, data=data[s.ix, , drop=FALSE], perturbations=perturbations[s.ix, , drop=FALSE], predn=predn)

		## performance estimation: R2
		mynet.r2 <- pred.score(data=data[s.ix, , drop=FALSE], pred=mynet.pred, method="r2")[nn]
		
		## performance estimation: NRMSE
		mynet.nrmse <- pred.score(data=data[s.ix, , drop=FALSE], pred=mynet.pred, method="nrmse")[nn]

		## performance estimation: MCC
		## discretize predicted gene expression data
		pred.discr <- data.discretize(data=mynet.pred, cuts=cuts.discr[dimnames(mynet.pred)[[2]]])
		mynet.mcc <- pred.score(data=mydata.discr[s.ix, dimnames(mynet.pred)[[2]], drop=FALSE], pred=pred.discr, method="mcc")[nn]

		## adjacency matrix
		topol <- regrnet2topo(net=mynet$net, coefficients=FALSE)

		## save results
		myr2 <- rbind(myr2, mynet.r2)
		mynrmse <- rbind(mynrmse, mynet.nrmse)
		mymcc <- rbind(mymcc, mynet.mcc)
		mytopo <- c(mytopo, list(topol))
	}
	dimnames(myr2)[[1]] <- dimnames(mynrmse)[[1]] <- dimnames(mymcc)[[1]] <- names(mytopo) <- names(mynets) <- paste("fold", 1:nfold, sep=".")
	
	## compute stability for each edge
	edgestab <- edgestab2 <- matrix(0, nrow=nrow(mytopo[[1]]), ncol=ncol(mytopo[[1]]), dimnames=dimnames(mytopo[[1]]))
	for(i in 1:length(mytopo)) {
		edgestab <- edgestab + mytopo[[i]]
	}
	edgestab <- edgestab / length(mytopo)
	## report stability of edges present in the global network
	edgestab2[mytopoglobal == 1] <- edgestab[mytopoglobal == 1]
	
	
	return(list("net"=mynetglobal, "net.cv"=mynets, "topology"=mytopoglobal, "topology.cv"=mytopo, "prediction.score.cv"=list("r2"=myr2, "nrmse"=mynrmse, "mcc"=mymcc), "edge.stability"=edgestab2, "edge.stability.cv"=edgestab))
}
