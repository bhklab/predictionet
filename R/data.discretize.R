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
	if(nrow(data) == 1) { res <- t(res) }
	dimnames(res) <- dimnames(data)
	return(res)
}
