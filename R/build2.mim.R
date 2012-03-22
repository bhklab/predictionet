### function to compute the matrix of mutual information (inspired from minet)
### returns a matrix containing the mutual information between all pairs of variables
## dataset: Matrix of continuous values (gene expressions for example); observations in rows, features in columns.
## estimator: type of correlation coefficient used to compute the mutual information
`.build2.mim` <- 
function(dataset, estimator=c("pearson", "spearman", "kendall")) {
	estimator <- match.arg(estimator)
	if( estimator=="pearson" || estimator=="spearman" || estimator=="kendall") {
		mim <- cor(dataset,method=estimator,use="complete.obs")^2
		diag(mim) <- 0
		maxi <- 0.999999
		mim[which(mim>maxi)] <- maxi
		mim <- -0.5*log(1-mim)
	} else { stop("unknown estimator") }
	
	mim[mim<0] <- 0
	return(mim)
}
