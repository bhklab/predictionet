### function computing the f1-score, comparing an inferred topology with a given topology
## net: inferred topology
## true.net: topology the user wants to compare the inferred topology with, e.g. the true network using generated datasets
### returns f1-score 
`eval.network` <- 
function(net, true.net) {
	
	tp<-fn<-fp<-fd<-0	
	
	### find counts for tp, fn, fp, fd
	for(i in 1:nrow(net)){
		for(j in 1:ncol(net)){
			if(true.net[i,j]>0){
				if(net[i,j]>0){
					tp<-tp+1
				}else{
					fn<-fn+1
				}
			}else if(true.net[j,i]>0 && net[i,j]>0){
				fd<-fd+1
			}else if( net[i,j]>0){
				fp<-fp+1
			}
		}
	}
	## considering edges in the wrong direction as false positives
	fp<-fp+fd
	
	## return f1 score
	return(2*tp/(2*tp+fn+fp))
}
