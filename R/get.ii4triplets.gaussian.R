### function to compute the interaction information of triplets of variables
### returns vector containing estimated interaction information values
## data: Matrix of continuous values (gene expressions for example); observations in rows, features in columns.
## triplets: matrix of triplets for which interaction information should be computed: rows containing triplets
`.get.ii4triplets.gaussian` <- 
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
