`.rank.genes.ensemble`<-function(models,data){
	
#	models<-res.ensemble$models
	mynames<-colnames(data)
	models.parents<-matrix(Inf,nc=ncol(models),nr=nrow(models))
	colnames(models.parents)<-models[1,]
	
	mim<- .build2.mim(data,estimator="pearson")
	for(i in 1:ncol(models)){
		for(j in 2:nrow(models)){
			models.parents[j,i] <- mim[models[j,i],models[1,i]]
			red<-0
			if(j>2){
				for(k in (j-1):2){
					red<-red+mim[models[k,i],models[j,i]]
				}
				red<-red/(j-2)
			}
			models.parents[j,i] <- models.parents[j,i] - red
		}
		
	}	
	return (list(score.local=models.parents))
}
