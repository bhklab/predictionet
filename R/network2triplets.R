### function finding all triplets in symmetric matrix X of the form A-B-C and A not connected with C
### X adjacency matrix (symmetric)
`.network2triplets` <- 
function(X) {
### step 1: find all pairs (matrix in which the first element has a smaller number than the second)
### step 2: identify the triplets using matrix of pairs
	if(sum(X)!=0){
		n <- dim(X)[1]
		m <- length(X[X!=0])
		
		step1 <- .get.pairs(X,n,m)
		
		pairs <- step1$pairs		 	# matrix of dimension m times 2, ordered pairs
		index <- step1$index		 	# index vector containing on position i the row number
		
		triplets <- matrix(0,nrow=90*m,ncol=3)
		cnt <- 1
		
		for(j in 1:m){
			
			p <- pairs[j,2]
			p2 <- pairs[j,1]
			
			p1 <- index[p]
			q1 <- index[p+1]-1
			
			indq1 <- 1
			while(q1==(-1)){
				q1 <- index[p+indq1]-1
				indq1 <- indq1+1
			}
			for(k in (p1:q1)){
				q2 <- pairs[k,2]
				
				p3 <- p2
				while(index[p3+1]==0){
					p3 <- p3+1
				}
				if(q2 > p2 && !((q2)%in% pairs[(index[p2]):(index[p3+1]-1),2])){ 
					triplets[cnt,] <- c(pairs[j,],q2)
					cnt <- cnt+1
				}
			}
		}
		return(triplets[triplets[,1]!=0,])
	}else{
		return(NULL)
	}		
}
