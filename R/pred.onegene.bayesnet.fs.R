### function predicting one gene at the time for the test data using the bayesian network model built on the training data
## net: bayesian network model
## data: Matrix of categorical values (discretized gene expressions for example); observations in rows, features in columns.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## predn: indices or names of variables (genes) to predict
### return caterigorical predictions for each gene
`.pred.onegene.bayesnet.fs` <- 
function(net, data, categories, perturbations, subset, predn) {
	#require(catnet)
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
