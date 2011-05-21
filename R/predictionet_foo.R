## functions of our netinf package

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
`netinf` <- 
function(data, categories, perturbations, priors, predn, priors.count=TRUE, priors.weight=0.5, maxparents=3, maxparents.force=FALSE, subset, method=c("regrnet", "regrnet.ensemble", "bayesnet", "bayesnet.ensemble"), regrmodel=c("linear", "linear.penalized"), causal=TRUE, seed=54321, retoptions="all", ...) {
	
	if(!missing(predn) && !is.null(predn) && (length(predn) < 2 && method!="regrnet.ensemble")) { stop("length of parameter 'predn' should be >= 2!") }
	maxparents <- ifelse(maxparents >= ncol(data), ncol(data)-1, maxparents)
	if(causal && maxparents > (ncol(data) * 0.5)) { warning("maximum number of parents may be too large, causal inference requires sparsity in the inferred network; please decrease maxparents parameter for better results!") }
	method <- match.arg(method)
	regrmodel <- match.arg(regrmodel)
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
		
		bnet <- fit.catnet(data=data, categories=categories, perturbations=perturbations, priors=priors, priors.weight=priors.weight, maxparents=maxparents, seed=seed, ...)
		
		if(retoptions=="all") {
		   return(list("method"=method, "net"=bnet))
		} else {
		   return(list("method"=method, "topology"=bayesnet2topo(net=bnet)))
		}
	}, 
	"bayesnet.ensemble"={
		stop("ensemble bayesian network inference is not implemented yet!")
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
		} else {
			if(max(priors, na.rm=TRUE) > 1 || min(priors, na.rm=TRUE) < 0) { stop("'priors' should contain probabilities of interactions if 'priors.count' is FALSE!") }
			## priors are probabilties taht should be rescale in [-1, 1]
			priors <- (priors - 0.5) * 2
		}
		bnet <- fit.regrnet.causal(data=data, perturbations=perturbations, priors=priors, predn=predn, maxparents=maxparents, priors.weight=priors.weight, causal=causal, regrmodel=regrmodel, seed=seed)
		if(retoptions=="all") {
		   return(list("method"=method, "net"=bnet))
		} else {
		   return(list("method"=method, "topology"=regrnet2topo(net=bnet,coefficients=FALSE), "topology.coeff"=regrnet2topo(net=bnet,coefficients=TRUE)))
		}
	}, 
	"regrnet.ensemble"={
		## fit an ensemble regression model
		rep_boot<-200
		maxnsol<-maxparents
		vec_ensemble <- .Call("mrmr_ensemble", data.matrix(data),maxparents, ncol(data), nrow(data), predn, length(predn), rep_boot, maxnsol, -1000)
		   net<-extract.adjacency.ensemble(data,vec_ensemble,predn)
		   models.equiv<-extract.all.parents(data,vec_ensemble,maxparents,predn)
		   return(list("method"=method, "net"=net, "models"=models.equiv))
	})
}

### Function predicting values of nodes from a previously inferred network
`netinf.predict` <- 
function(net, data, categories, perturbations, subset, predn) {
	res <- matrix(NA, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	switch(net$method, 
	"bayesnet"={
		res <- pred.onegene.bayesnet.fs(net=net$net, data=data, categories=categories, perturbations=perturbations, subset=subset, predn=predn)
	}, 
	"regrnet"={
		res <- pred.onegene.regrnet.fs(net=net$net, data=data, perturbations=perturbations, subset=subset, predn=predn)
	}, 
	"regrnet.ensemble"={
		## fit an ensemble regression model
		stop("Ensemble regression-based network inference is not implemented yet!")
	})
	return(res)
}

### function fitting a bayesian network model
## data: Matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## priors: matrix of prior information avalable for gene-gene interaction (parents in columns, children in rows). Values may be probabilities or any other weights (citations count for instance)
## maxparents: maximum number of parents allowed for each gene
## seed: seed to get deterministic results
### returns a bayesian network model
'fit.catnet' <- 
function(data, categories, perturbations, priors, priors.weight, maxparents=3, maxparents.force=FALSE, seed=54321, ...) {
	## run catnet
	require(catnet)
	catnet::cnSetSeed(seed)
	## be aware that catnet package consider any adjacency matrix to have parents in COLUMNS and children in ROWS, that is the inverse of the predictionet package
	## topology identified from the priors
	priorparents <- sapply(1:ncol(priors), function(x, y) { tt <- which(y[ , x] > 0.5); if(length(tt) < 1) { tt <- NULL }; return(list(tt)); }, y=adj.remove.cycles(adjmat=priors)[[1]])
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
		return(list("varnames"=colnames(data), "input"=myparents, "model"=ee))
	} else {
		## use data and priors infer the best topology
		if(priors.weight == 0) {
			## do not use the priors for identifying the network topology
			## uniform priors
			priors[] <- 0.5
		}
		## seed the network inference with the prior topology
		#priororder <- catnet::cnOrder(priorparents)
		#ee.prior <- catnet::cnSearchOrder(data=t(data), perturbations=t(perturbations), maxParentSet=maxparents, nodeOrder=priororder, edgeProb=t(priors), ...)
		ee.prior <- NULL
		ee <- catnet::cnSearchSA(data=t(data), nodeCats=categories, perturbations=t(perturbations), selectMode="BIC", maxParentSet=maxparents, priorSearch=ee.prior, edgeProb=t(priors), dirProb=t(priors), echo=FALSE, ...)
		if(maxparents.force) { ee <- ee@nets[order(ee@complexity, decreasing=TRUE)[1]] } else { ee <- cnFindBIC(ee) }
		myparents <- lapply(ee@parents, function(x, y) { if(!is.null(x)) { x <- y[x]}; return(x); }, y=ee@nodes)
		names(myparents) <- ee@nodes
		return(list("varnames"=colnames(data), "input"=myparents, "model"=ee))
	}
}

### function fitting a regression model for each gene in the data
## data: Matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## priors: matrix of prior information avalable for gene-gene interaction (parents in columns, children in rows). Values may be probabilities or any other weights (citations count for instance)
## predn: indices or names of genes to fit during network inference. If missing, all the genes will be used for network inference
## priors.weight: real value in [0, 1] specifying the weight to put on the priors (0=only the data are used, 1=only the priors are used to infer the topology of the network). If 'priors.weight' is missing it will be optimized gene by hene in an automatic way.
## maxparents: maximum number of parents allowed for each gene
## maxparents.force: if TRUE it forces the method to select the maximum number of parents for each variable in the network, FALSE otherwise. 
### returns a regression network 
`fit.regrnet.causal` <- 
function(data, perturbations, priors, predn, maxparents=3, maxparents.force=FALSE, priors.weight=0.5, regrmodel=c("linear", "linear.penalized"), causal=TRUE, seed=54321) {
	
	if(!missing(seed)) { set.seed(seed) }
	if(missing(predn) || is.null(predn) || length(predn) == 0) { predn <- 1:ncol(data) } else { if(is.character(predn)) { predn <- match(predn, dimnames(data)[[2]]) } else { if(!is.numeric(predn) || !all(predn %in% 1:ncol(data))) { stop("parameter 'predn' should contain either the names or the indices of the variables to predict!")} } }
	nn <- dimnames(data)[[2]][predn] ## variables (genes) to fit during network inference
	regrmodel <- match.arg(regrmodel)
	regrnet <- NULL
	
	########################
	### data preparation
	########################
	dd <- data.frame(data)
	diag(priors) <- 0 ## no self-loops
	
	########################
	### find ranked parents for all genes
	########################
	mat.ranked.parents <- rank.genes.causal.perturbations(priors=priors, data=dd, perturbations=perturbations, predn=predn, priors.weight=priors.weight, maxparents=maxparents, maxparents.force, causal=causal)
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
			switch(regrmodel, 
			"linear"={
				mm.reg <- lm(formula=ff, data=dd)
			},
			"linear.penalized"={
				require(penalized)
				optlambda1 <- optL1(response=ff, data=dd, model="linear", lambda2=0, minlambda1=1, maxlambda1=10, trace=FALSE, fold=10)
				mm.reg <- penalized(response=ff, data=dd, model="linear", lambda1=optlambda1$lambda, lambda2=0, trace=FALSE)
				## cor(data[ , nn[i]], predict(object=mm.reg, data=dd)[ , 1])
				## compare with lm
				## mm <- lm(formula=ff, data=dd)
				## cor(data[ , nn[i]], predict(mm, newdata=dd))
			})
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
function(priors, data, perturbations, predn, priors.weight, maxparents, maxparents.force=FALSE, causal=TRUE) {

	########################
	### initialize variables
	########################
	res <- matrix(NA, nrow=length(predn), ncol=ncol(data) - 1, dimnames=list(dimnames(data)[[2]][predn], seq(1, ncol(data) - 1)))
	res.netw <- matrix(0, nrow=ncol(data), ncol=length(predn), dimnames=list(dimnames(data)[[2]], dimnames(data)[[2]][predn]))
	
	data.original <- data.matrix(data)
	
	if(priors.weight != 1) {
		## data will be used to rank the potential parent genes
		mrmr.global <- NULL
		if(all(apply(perturbations, 2, sum, na.rm=TRUE) == 0) & sum(is.na(perturbations)) == 0) {
			## no perturbations
			## compute mrmr score for the target genes
			## mrmr_adapted(SEXP data, SEXP maxparents, SEXP nvar, SEXP nsample, SEXP predn, SEXP npredn, SEXP threshold);
			mrmr.global <- .Call("mrnet_adapted", data.original, maxparents, ncol(data.original), nrow(data.original), predn, length(predn), -1000)
			mrmr.global <- matrix(mrmr.global, nrow=ncol(data.original), ncol=ncol(data.original), byrow=FALSE)
			dimnames(mrmr.global) <- list(colnames(data.original), colnames(data.original))
			mrmr.global[mrmr.global == -1000] <- NA
		}

		for(i in 1:length(predn)) {
			data <- data.original
			## remove observation where the target variabl ewas perturbed (no link between observation and target variable)
			ind <- which(is.na(perturbations[, predn[i]]) | perturbations[, predn[i]] == 1)
			if(length(ind)>0) {
				data <- data[-ind,]
				## since the target gene has been perturbed in some of the experiments we should recompute the mim and mrmr matrices
				## compute mrmr score for the target genes
				## mrmr_adapted(SEXP data, SEXP maxparents, SEXP nvar, SEXP nsample, SEXP predn, SEXP npredn, SEXP threshold);
				mrmr <- .Call("mrnet_adapted", data, maxparents, ncol(data), nrow(data), predn, length(predn), -1000)
				mrmr <- matrix(mrmr, nrow=ncol(data), ncol=ncol(data), byrow=FALSE)
				dimnames(mrmr) <- list(colnames(data),colnames(data))
				mrmr[mrmr == -1000] <- NA
			} else { mrmr <- mrmr.global }
			## mrmr scores on the diagonal should be NA
			diag(mrmr) <- NA
			## select only the variables specified by the user (predn)
			mrmr <- mrmr[ , predn, drop=FALSE]

			########################
			### compute mrmr score matrix
			########################
			## normalization of mrmr scores
			mrmr <- mrmr/max(abs(mrmr), na.rm=TRUE)
			## avoid mrmr score of exactly 0
			mrmr[!is.na(mrmr) & mrmr == 0] <- 1e-16
			score <- matrix(0, nrow=ncol(data), ncol=ncol(data), dimnames=list(colnames(data), colnames(data)))
			score[ , predn][!is.na(mrmr)] <- mrmr[!is.na(mrmr)]

			if(causal) {
				########################
				## find all adjacencies of node to predict in score
				########################
				tmp.adj <- matrix(0,ncol=ncol(data),nrow=ncol(data), dimnames=list(colnames(data), colnames(data)))
				tmp.adj[predn[i],] <- tmp.adj[ ,predn[i]] <- score[ ,predn[i]]

				########################
				## determine all triplets X.k-X.j-X.l from these adjecencies
				########################
				trip <- (network2triplets(tmp.adj))
				parents <- NULL
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
					if(nrow(trip)>0){
						for(j in 1:nrow(trip)){
							if(!is.nan(trip[j,4]) && trip[j,4]<0) {
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
			## parents into network matrix
			if(!is.null(parents)) {
				## take the minimum between the mrmr scores and 0 for the genes for which we were not able to infer causality
				tt <- score[ ,predn[i]]
				parix <- which(score[ ,predn[i]] != 0)
				tt[setdiff(parix, parents)] <- ifelse(tt[setdiff(parix, parents)] < 0, tt[setdiff(parix, parents)], 0)
				## it is not smooth but it penalized the genes with positive mrmr score but for which we do not evidence for causality
				res.netw[ ,i] <- tt
			}
		}
	}
	
	########################
	## (1-w)*M+w*P weighted score from prior knowledge and data
	########################
	score <- (1-priors.weight)*res.netw+priors.weight*priors[ , predn, drop=FALSE]
	## no self-loops
	diag(score[dimnames(score)[[2]], dimnames(score)[[2]]]) <- 0
	for(j in 1:length(predn)){
		tmp.s <- sort(score[,j],decreasing=TRUE,index.return=TRUE)
		if(priors.weight != 1 && maxparents.force) { ind.rm <- order(tmp.s$x, decreasing=TRUE)[1:maxparents] } else { ind.rm <- (which(tmp.s$x<=0)) }
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
	if(missing(predn) || is.null(predn) || length(predn) == 0) { predn <- match(names(net$input), dimnames(data)[[2]]) } else { if(is.character(predn)) { predn <- match(predn, dimnames(data)[[2]]) } else { if(!is.numeric(predn) || !all(predn %in% 1:ncol(data))) { stop("parameter 'predn' should contain either the names or the indices of the variables to predict!")} } }
	if(!all(is.element(predn, match(names(net$input), dimnames(data)[[2]])))) { stop("some genes cannot be predicted because they have not been fitted in the network!")}
	## variables to predict
	nnix <- dimnames(data)[[2]][predn]
	
	## matrix to store the predictions
	preds <- matrix(NA, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	for (i in 1:length(nnix)) {
		## extract the fitted model
		model.i <- net$model[[nnix[i]]]
		if(!is.numeric(model.i)) {
			## non null model
			switch(class(model.i),
			"lm"={
				## linear model
				preds[ , nnix[i]] <- predict(object=model.i, newdata=data.frame(data))
			},
			"penfit"= {
				## penalized linear model
				require(penalized)
				preds[ , nnix[i]] <- predict(object=model.i, data=data.frame(data))[ , 1]
			})
		} else { preds[ , nnix[i]] <- model.i }
	}
	preds[perturbations] <- NA
	return(preds)
}

### function predicting one gene at the time for the test data using the bayesian network model built on the training data
## net: bayesian network model
## data: Matrix of categorical values (discretized gene expressions for example); observations in rows, features in columns.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## predn: indices or names of variables (genes) to predict
### return caterigorical predictions for each gene
`pred.onegene.bayesnet.fs` <- 
function(net, data, categories, perturbations, subset, predn) {
	require(catnet)
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
	if(missing(predn) || is.null(predn) || length(predn) == 0) { predn <- match(names(net$input), dimnames(data)[[2]]) } else { if(is.character(predn)) { predn <- match(predn, dimnames(data)[[2]]) } else { if(!is.numeric(predn) || !all(predn %in% 1:ncol(data))) { stop("parameter 'predn' should contain either the names or the indices of the variables to predict!")} } }
	if(!all(is.element(predn, match(names(net$input), dimnames(data)[[2]])))) { stop("some genes cannot be predicted because they have not been fitted in the network!")}	
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
	## variables to predict
	nnix <- dimnames(data)[[2]][predn]
	## matrix to store the predictions
	preds <- matrix(NA, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	for (i in 1:length(nnix)) {
		datat <- data
		datat[ , nnix[i]] <- NA
		## predict the missing values
		preds[ , nnix[i]] <- t(catnet::cnPredict(object=net$model, data=t(datat)))[ , nnix[i]]
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
## net: bayesian network
### returns adjacency matrix: parents in rows and children in columns (res[i,j]==1 means edge from i to j)
`bayesnet2topo` <- 
function(net) {
	require(catnet)
	return(t(catnet::cnMatParents(net$model)))
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
		res <- bayesnet2topo(net=net$net)
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
function(data, categories, perturbations, priors, predn, priors.count=TRUE, priors.weight=0.5, maxparents=3, maxparents.force=FALSE, subset, method=c("regrnet", "regrnet.ensemble", "bayesnet", "bayesnet.ensemble"), regrmodel=c("linear", "linear.penalized"), nfold=10, causal=TRUE, seed, retoptions="all", ...) {
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
		mydata.discr <- data.discretize(data=data, cuts=cuts.discr)
		categories <- cuts.discr
	} else { mydata.discr <- data }
	
	## infer network from the whole dataset
	mynetglobal <- netinf(data=data, categories=categories, perturbations=perturbations, priors=priors, predn=predn, priors.count=priors.count, priors.weight=priors.weight, maxparents=maxparents, method=method, regrmodel=regrmodel, causal=causal)
	## and the global topology
	mytopoglobal <- net2topo(net=mynetglobal)
	
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
		mynet <- netinf(data=data[-s.ix, , drop=FALSE], categories=categories, perturbations=perturbations[-s.ix, , drop=FALSE], priors=priors, predn=predn, priors.count=priors.count, priors.weight=priors.weight, maxparents=maxparents, method=method, regrmodel=regrmodel, causal=causal)
		mynets <- c(mynets, list(mynet))
		
		## compute predictions
		mynet.pred <- netinf.predict(net=mynet, data=data[s.ix, , drop=FALSE], categories=categories, perturbations=perturbations[s.ix, , drop=FALSE], predn=predn)
		
		resnull <- rep(NA, ncol(data))
		names(resnull) <- colnames(data)
		## performance estimation: R2
		mynet.r2 <- resnull
		if(method %in% c("regrnet")) { mynet.r2 <- pred.score(data=data[s.ix, , drop=FALSE], pred=mynet.pred, method="r2")[nn] }
		## performance estimation: NRMSE
		mynet.nrmse <- resnull
		if(method %in% c("regrnet")) { mynet.nrmse <- pred.score(data=data[s.ix, , drop=FALSE], pred=mynet.pred, method="nrmse")[nn] }
		## performance estimation: MCC
		pred.discr <- mynet.pred
		if(method %in% c("regrnet")) {
			## discretize predicted gene expression data 
			pred.discr <- data.discretize(data=mynet.pred, cuts=cuts.discr[dimnames(mynet.pred)[[2]]])	
		}
		mynet.mcc <- pred.score(data=mydata.discr[s.ix, dimnames(mynet.pred)[[2]], drop=FALSE], pred=pred.discr, method="mcc")[nn]

		## adjacency matrix
		topol <- net2topo(net=mynet)

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
	
	if(retoptions=="all") {
		return(list("net"=mynetglobal, "net.cv"=mynets, "topology"=mytopoglobal, "topology.cv"=mytopo, "prediction.score.cv"=list("r2"=myr2, "nrmse"=mynrmse, "mcc"=mymcc), "edge.stability"=edgestab2, "edge.stability.cv"=edgestab))
	} else {
		mytopo2<-NULL
		for(i in 1:nfold){
			mytopo2<-c(mytopo2, list(net2topo(net=mynets[[i]], coefficients=TRUE)))
		}
		names(mytopo2)<- paste("fold", 1:nfold, sep=".")
		return(list("topology.coeff"=net2topo(net=mynetglobal, coefficients=TRUE), "topology.cv.coeff"=mytopo2, "topology"=mytopoglobal, "topology.cv"=mytopo, "prediction.score.cv"=list("r2"=myr2, "nrmse"=mynrmse, "mcc"=mymcc), "edge.stability"=edgestab2, "edge.stability.cv"=edgestab))
	}
}

'extract.all.parents' <-
function(data,res.main,maxparents,predn) {
### function taking the output of the regrnet.ensemble method and returns a matrix
### containing one equivalent model in each column. The target variable is in the first row.
### res.main: 	output of regrnet.ensemble
### maxparents:	maxparents parameter of the netinf method
### predn:	list of target variables for which ensemble method was run.
	
	final<-NULL
	cnt_main<-1
	for(imain in 1:length(predn)){
		res.vec<-res.main[cnt_main:(cnt_main+2*res.main[cnt_main])]
		if(length(res.vec)>3){
		cnt_main<-cnt_main+2*res.main[cnt_main]+1
		nsol<-sum(res.vec==0)
		res<-matrix(0,nc=nsol,nr=(maxparents+1))
		
		val<-res.vec[2:(res.vec[1]+1)]
		ind<-res.vec[(res.vec[1]+2):(2*res.vec[1]+1)]
		res[1,1:ind[1]]<-rep(val[1],ind[1])
		nvar<-length(val)
		level<-2
		cnt<-1
		nelem<-0
		sum_old<-sum(ind[1])
		last_level<-FALSE
		i<-2
		ind2<-1
		
		while(i<=nvar && !last_level){
			if(ind[i]!=0){
				if(ind[i]>1){
					tmp<-res[,(ind2+1):nsol]
					res[level,ind2]<-val[i]
					for(j in 1:level){
						res[j,(ind2+1):(ind2+ind[i]-1)]<-rep(res[j,ind2],(ind[i]-1))
					}
					if((nsol-(ind2+ind[i]-1))>0){
						res<-cbind(res[,1:(ind2+ind[i]-1)],tmp[,1:(nsol-(ind2+ind[i]-1))])
					}else{
						res<-res[,1:(ind2+ind[i]-1)]
					}
				}
				res[level,ind2:(ind2-1+ind[i])]<-rep(val[i],ind[i])
				ind2<-ind2+ind[i]
			}else{
				res[level,ind2]<-val[i]
				ind2<-ind2+1
			}
			nelem<-nelem+ind[i]
			if(cnt==sum_old){			
				sum_old<-nelem
				if(nelem==0){
					last_level<-TRUE
				}
				nelem<-0
				cnt<-1
				level<-level+1
				ind2<-1
			}
			else if(cnt< sum_old){
				cnt<-cnt+1
			}
			i<-i+1
		}
		final<-cbind(final,res)
		}
	}
	dimension<-dim(final)
	final<-colnames(data)[final]
	dim(final)<-dimension
	return(final)
}

'extract.adjacency.ensemble' <- 
function(data,res.vec,predn) {
	res<-matrix(0,nc=ncol(data),nr=ncol(data),dimnames=list(colnames(data),colnames(data)))
	ind.start<-1
	for(i in 1:length(predn)){
		target<-res.vec[ind.start+1]
		for(j in ((ind.start+2):(ind.start+res.vec[ind.start]))){
			res[res.vec[j],target]<-1
			
		}
		ind.start<-ind.start+2*res.vec[ind.start]+1
	}
	return(res)
}

## function to remove cycles in an adjacency matrix
## adjmat: adjacency matrix; parents in rows, children in collumns
## return the adjacency matrix with the lowest entries being removed such that the network is now acyclic
'adj.remove.cycles' <- 
function(adjmat) {
	
	####################
	## internal functions
	####################
	## adjmat: adjacency matrix
	## nodest: status of nodes (vector); 0 if not visited, 1 otherwise
	## nodix.from: index of the previous node
	## nodix: index of the current node
	'.adj.remove.cycles.DFS' <- 
	function(adjmat, adjmat.mask, nodix.from, nodix, nodest) {

		if(nodix.from %in% (1:length(nodest)) && nodest[nodix] == 1) { ## the current node has already been visited
			## remove edge
			adjmat.mask[nodix.from, nodix] <- TRUE
			return(adjmat.mask)
		} else {
			## the current node has never been visited
			## remove the edges previously identified as responsible for a cycle
			adjmat2 <- adjmat
			adjmat2[adjmat.mask] <- 0
			ee <- which(adjmat2[nodix, ] > 0) ## all children of nodix
			ee <- ee[order(adjmat2[nodix, ee], decreasing=FALSE)] ## order them based on number of citations
			if(length(ee) > 0) { ## some children to look at
				for(i in 1:length(ee)) { ## for all children
					nodest2 <- nodest
					nodest2[nodix] <- 1
					adjmat.mask <- .adj.remove.cycles.DFS(adjmat=adjmat, adjmat.mask=adjmat.mask, nodix.from=nodix, nodix=ee[i], nodest=nodest2)
				}
				return(adjmat.mask)
			} else { return(adjmat.mask) }
		}
	}
	####################
	if(class(adjmat) != "matrix" && nrow(adjmat) != ncol(adjmat)) { stop("the adjacency matrix should be square!") }
	nNames <- dimnames(adjmat)[[1]]
	adjmat <- adjmat[nNames, nNames, drop=FALSE]
	adjmat2 <- adjmat
	adjmat.mask <- matrix(FALSE, nrow=nrow(adjmat), ncol=ncol(adjmat), dimnames=dimnames(adjmat))
	
	## remove bidirectional edges wrt the number of interactions
	## on the diagonal
	iix <- which(diag(adjmat2) != 0)
	if(length(iix) > 0) {
		## some entries in the diagonal should be removed
		adjmat.mask[cbind(iix, iix)] <- TRUE
	}
	adjmat2[adjmat.mask] <- 0
	## bidirectional edges
	iix <- which(upper.tri(adjmat2), arr.ind=TRUE)
	iix2 <- t(apply(iix, 1, function(x, y) { if(y[x[1], x[2]] <= y[x[2], x[1]]) { return(x) } else { return(rev(x)) } }, y=adjmat2))
	iix2 <- iix2[which(adjmat2[iix2] != 0), ]
	adjmat.mask[iix2] <- TRUE
	adjmat2[adjmat.mask] <- 0
	
	## order nodes with the maximum prior
	nnix <- order(apply(adjmat2, 1, max), decreasing=FALSE)
	for(i in 1:length(nnix)) {
		nodest <- rep(0, ncol(adjmat2))
		names(nodest) <- colnames(adjmat2)
		adjmat.mask <- .adj.remove.cycles.DFS(adjmat=adjmat2, adjmat.mask=adjmat.mask, nodix.from=-1, nodix=nnix[i], nodest)
		adjmat2[adjmat.mask] <- 0
	}
	
	return(list("adjmat.acyclic"=adjmat2, "adjmat.removed"=adjmat.mask))
}
