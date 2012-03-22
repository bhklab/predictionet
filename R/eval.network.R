### function computing the f1-score, comparing an inferred topology with a given topology
## topo: inferred topology
## true.topo: topology the user wants to compare the inferred topology with, e.g. the true topowork using generated datasets
### returns f1-score 
`eval.network` <- 
function(topo, true.topo) {
	
	tp <- fn <- fp <- fd <- 0	
	
	### find counts for tp, fn, fp, fd
	for(i in 1:nrow(topo)){
		for(j in 1:ncol(topo)){
			if(true.topo[i,j]>0){
				if(topo[i,j]>0){
					tp <- tp+1
				}else{
					fn <- fn+1
				}
			}else if(true.topo[j,i]>0 && topo[i,j]>0){
				fd <- fd+1
			}else if(topo[i,j]>0){
				fp <- fp+1
			}
		}
	}
	## considering edges in the wrong direction as false positives
	fp <- fp+fd
	
	## return f1 score
	return(2*tp/(2*tp+fn+fp))
}
