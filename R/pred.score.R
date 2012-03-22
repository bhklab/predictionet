### function computing perfomance of prediction; methods include rnmse and mcc
### user has to provide data/pred accordingly (either both continuous or both discretized)
## data: matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## pred: Matrix of continuous or categorical values: predictions in rows, features in columns
## categories: list of categories for each of the nodes (genes) in 'data' matrix. Categories should be integers from 1 to n.
## method: scoring metric: nrmse, mcc,...
`pred.score` <- 
function(data, pred, categories, method=c("r2", "nrmse", "mcc")) {
	if(is.vector(data)) { data <- matrix(data, ncol=1, dimnames=list(names(data), "X")) }
	if(is.vector(pred)) { pred <- matrix(pred, ncol=1, dimnames=list(names(pred), "X")) }
	method <- match.arg(method)
	data.new <- matrix(0,ncol=ncol(pred), nrow=nrow(data), dimnames=list(rownames(data), colnames(pred)))
	for(i in 1:ncol(pred)){
		data.new[,i] <- data[,colnames(pred)[i]]
	}
	data <- data.new
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
			   })
		return(myperf)
	}
	
	if(is.vector(data) && is.vector(pred)) {
		if(length(data) != length(pred)) { stop("length of 'data' and 'pred' should be equal!") }
		myperf <- myfoo(data=data, pred=pred, nbcat=nbcat, method=method)
	} else {
		if(!all(dim(data) == dim(pred))) { stop("dimensions of 'data' and 'pred' should be equal!") }
		myperf <- apply(X=rbind(nbcat, data, pred), MARGIN=2, FUN=function(x, foo) {
			nbcat <- x[1]
			ll <- (length(x)-1)/2
			data <- x[2:(ll+1)]
			pred <- x[(ll+2):((2*ll)+1)]
			return(foo(data=data, pred=pred, nbcat=nbcat, method=method))
		}, foo=myfoo)
	}
	
	return(myperf)
}
