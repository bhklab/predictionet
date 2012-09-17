### function finding all triplets in symmetric matrix X of the form A-B-C and A not connected with C
### X adjacency matrix (symmetric)
`.network2triplets` <- 
function(X,target) {
### step 1: find all pairs (matrix in which the first element has a smaller number than the second)
### step 2: identify the triplets using matrix of pairs
	if(sum(X)!=0){
		cand_dim1<-which(X[,target]>0)
		cand_dim2<-which(X[target,]>0)
		cand<-union(cand_dim1,cand_dim2)
		cand<-sort(cand)
		ntrip<-(length(cand)-1)*length(cand)/2
		if(ntrip>0){
			trip<-matrix(0,ncol=3,nrow=ntrip)
			trip[,2]<-rep(target,ntrip)
			cnt<-1
			to.rm<-NULL
			for(i in 1:(length(cand)-1)){
				for(j in (i+1):(length(cand))){
					trip[cnt,1]<-cand[i]
					trip[cnt,3]<-cand[j]
					if(X[cand[i],cand[j]]>0){
						to.rm<-c(to.rm,cnt)
					}
					cnt<-cnt+1
					
				}
			}
			
			if(!is.null(to.rm)){
				if(length(to.rm)==ntrip){
					return(NULL)
				}else{
					return(trip[-to.rm,,drop=FALSE])
				}
			}else{
				return(trip)
			}
		}else{
			return(NULL)
		}
		
	}
}
