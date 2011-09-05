### function predicting one gene at the time for the test data using the regression model built on the training data
## net: regression model for each gene
## data: Matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## predn: indices or names of variables (genes) to predict
### return continuous predictions for each gene
`.pred.onegene.regrnet.fs` <- 
function(net, data, perturbations, subset, predn) {
	if(missing(perturbations) || is.null(perturbations)) {
		perturbations <- matrix(FALSE, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	} else {
		if(nrow(perturbations) == 1) {
			perturbations[1, ] <- as.logical(perturbations[1, ])
		} else { perturbations <- apply(perturbations, 2, as.logical) }
		dimnames(perturbations) <- dimnames(data)
	}
	
	if(!missing(subset)) {
		data <- data[subset, , drop=FALSE]
		perturbations <- perturbations[subset, , drop=FALSE]
	}
	if(missing(predn) || is.null(predn) || length(predn) == 0) { predn <- match(colnames(net), dimnames(data)[[2]]) } else { if(is.character(predn)) { predn <- match(predn, dimnames(data)[[2]]) } else { if(!is.numeric(predn) || !all(predn %in% 1:ncol(data))) { stop("parameter 'predn' should contain either the names or the indices of the variables to predict!")} } }
	if(!all(is.element(predn, match(colnames(net), dimnames(data)[[2]])))) { stop("some genes cannot be predicted because they have not been fitted in the network!")}
## variables to predict
	nnix <- dimnames(data)[[2]][predn]
	
## matrix to store the predictions
	preds <- matrix(NA, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	preds[,colnames(net )]<-t(t(net)%*% t(cbind(rep(1,nrow(data)),data)))
	
	preds[perturbations] <- NA
	return(preds)
}