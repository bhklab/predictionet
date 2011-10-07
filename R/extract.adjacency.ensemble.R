'.extract.adjacency.ensemble' <- 
function(data,res.vec,predn) {
	res <- matrix(0,nc=ncol(data),nr=ncol(data),dimnames=list(colnames(data),colnames(data)))
	ind.start <- 1
	for(i in 1:length(predn)){
		target <- res.vec[ind.start+1]
		for(j in ((ind.start+2):(ind.start+res.vec[ind.start]))){
			res[res.vec[j],target] <- 1
			
		}
		ind.start <- ind.start+2*res.vec[ind.start]+1
	}
	return(res)
}
