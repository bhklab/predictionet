### function fitting a regression model for each gene in the data
## data: Matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## priors: matrix of prior information avalable for gene-gene interaction (parents in columns, children in rows). Values may be probabilities or any other weights (citations count for instance)
## predn: indices or names of genes to fit during network inference. If missing, all the genes will be used for network inference
## priors.weight: real value in [0, 1] specifying the weight to put on the priors (0=only the data are used, 1=only the priors are used to infer the topology of the network). If 'priors.weight' is missing it will be optimized gene by hene in an automatic way.
## maxparents: maximum number of parents allowed for each gene
## maxparents.push: if TRUE it forces the method to select the maximum number of parents for each variable in the network, FALSE otherwise. 
### returns a regression network 
`.fit.regrnet.causal` <- 
function(data, perturbations, priors, predn, maxparents=3, maxparents.push=FALSE, priors.weight=0.5, regrmodel=c("linear", "linear.penalized"), causal=TRUE, seed=54321) {
	if(causal && maxparents > (ncol(data) * 0.5)) { warning("maximum number of parents may be too large, causal inference requires sparsity in the inferred network; please decrease maxparents parameter for better results!") }
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
	rrnetw <- .rank.genes.causal.perturbations(priors=priors, data=dd, perturbations=perturbations, predn=predn, priors.weight=priors.weight, maxparents=maxparents, causal=causal, maxparents.push=maxparents.push)
	mat.ranked.parents <- rrnetw$parents
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
				   #require(penalized)
				   optlambda1 <- penalized::optL1(response=ff, data=dd, model="linear", lambda2=0, minlambda1=1, maxlambda1=10, trace=FALSE, fold=10)
				   mm.reg <- penalized::penalized(response=ff, data=dd, model="linear", lambda1=optlambda1$lambda, lambda2=0, trace=FALSE)
					## cor(data[ , nn[i]], predict(object=mm.reg, data=dd)[ , 1])
					## compare with lm
					## mm <- lm(formula=ff, data=dd)
					## cor(data[ , nn[i]], predict(mm, newdata=dd))
				   })
		}
		regrnet <- c(regrnet, list(mm.reg))
	}
	names(regrnet) <- nn
	
	return(list("varnames"=colnames(data), "input"=myparents, "model"=regrnet, "edge.relevance"=rrnetw$score))
}
