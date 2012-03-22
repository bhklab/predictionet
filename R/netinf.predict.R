### Function predicting values of nodes from a previously inferred network
## net: object returned by netinf 
## data: matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## categories: Either an integer specifying the number of categories used to discretize the variables or a list of categories for each of the variables in 'data' matrix. Categories should be integers from 1 to n. If method='bayesnet', this parameter should be specified by the user.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## subset: vector of indices to select only subset of the observations
## predn: indices or names of genes to fit during network inference. If missing, all the genes will be used for network inference
## returns matrix containg a prediction for each sample and each variable in predn

`netinf.predict` <- 
function(net, data, categories, perturbations, subset, predn, method=c("linear", "linear.penalized", "cpt")) {
	method <- match.arg(method)
	if(length(method) > 1) { stop("only one prediction can be specified!") }
	if(is.element(method, c("linear", "linear.penalized")) && !is.element("lrm", names(net))) { stop(sprintf("the network model should be made prediction first!\nuse function net2pred with method %s", method)) }
	if(is.element(method, c("cpt")) && !is.element("cpt", names(net))) { stop(sprintf("the network model should be made prediction first!\nuse function net2pred with method %s", method)) }
	if(missing(perturbations) || is.null(perturbations)) {
		perturbations <- matrix(FALSE, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	} else {
		if(nrow(perturbations) == 1) {
			perturbations[1, ] <- as.logical(perturbations[1, ])
		} else { perturbations <- apply(perturbations, 2, as.logical) }
		dimnames(perturbations) <- dimnames(data)
	}
	ensemble <- net$ensemble
	res <- matrix(NA, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	switch(method, 
	"cpt"={
			stop("cpt method is not implemented yet!")
			#res <- .pred.onegene.bayesnet.fs(net=net$net, data=data, categories=categories, perturbations=perturbations, subset=subset, predn=predn)
	}, 
	"linear"=, 
	"linear.penalized"={
		   if(ensemble) {
				if(missing(predn)) { predn <- colnames(net$topology) } else if(length(intersect(predn,colnames(net$topology)))>0) { predn <- match(predn,colnames(data)) }
		   } else {
				if(missing(predn)) { predn <- seq(1,ncol(data)) } else if(length(intersect(predn,colnames(data)))>0) { predn <- match(predn,colnames(data)) }
		   }
			## make the predictions
			
#	res <- .pred.onegene.regrnet.fs(topo.coeff=.regrnet2matrixtopo(net), data=data, perturbations=perturbations, subset=subset, predn=predn, ensemble=ensemble)
		   res <- .pred.onegene.regrnet.fs(topo.coeff=net$lrm, data=data, perturbations=perturbations, subset=subset, predn=predn, ensemble=ensemble)

	}, 
	stop("no default parameter for method!")
	)
	return(res)
}
