### function returning an ordered matrix containing all pairs A-B in the symmetrix matrix X, 
### returns vector, element i contains position of start node i-... in matrix of pairs
### input, symmetrix matrix X containing 0 and 1, number of rows/columns, number of pairs 
### (=number of elements equal to 1)
`.get.pairs` <- 
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
