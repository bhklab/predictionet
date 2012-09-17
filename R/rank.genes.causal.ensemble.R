`.rank.genes.causal.ensemble`<-function(models,data){
	
#	models<-res.ensemble$models
	mynames<-colnames(data)
	models.parents<-matrix(Inf,ncol=ncol(models),nrow=nrow(models))
	colnames(models.parents)<-models[1,]
	score.causal <- matrix(Inf, nrow=ncol(data), ncol=ncol(data), dimnames=list(colnames(data), colnames(data)))
	
	for(i in 1:ncol(models)){
		tmp.adj <- matrix(0,ncol=ncol(data),nrow=ncol(data), dimnames=list(colnames(data), colnames(data)))
		tmp.adj[models[1,i],models[2:nrow(models),i]]<-tmp.adj[models[2:nrow(models),i],models[1,i]]<-1
		
		trip <- (.network2triplets(tmp.adj,models[1,i]))
		
		if(length(trip)>0){
			if(length(trip)==3){
				trip <- (t(as.matrix(trip)))
			}
########################
## compute I(X.k,X.l)-I(X.k,X.l|X.j)
########################
			
			trip <- as.data.frame(cbind(trip,(.get.ii4triplets.gaussian(data,trip))))
			
			if(nrow(trip)>0){
				for(j in 1:nrow(trip)){
					
					ind1<-mynames[trip[j,1]]
					ind3<-mynames[trip[j,3]]
					
					ind1<-which(models[,i]==ind1)
					ind3<-which(models[,i]==ind3)
					
					if(!is.nan(trip[j,4])){
						if( trip[j,4] < models.parents[ind1,i]){
							models.parents[ind1,i]<-trip[j,4]
						}
						if( trip[j,4] < models.parents[ind3,i]){
							models.parents[ind3,i]<-trip[j,4]
						}
					}
					
					if(!is.nan(trip[j,4]) && trip[j,4] < score.causal[trip[j,1],trip[j,2]]){
						score.causal[trip[j,1],trip[j,2]] <- trip[j,4]
					}
					if(!is.nan(trip[j,4]) && trip[j,4] < score.causal[trip[j,3],trip[j,2]]){
						score.causal[trip[j,3],trip[j,2]] <- trip[j,4]
					}
				}
			}
		}
	}	
	models.parents[2:nrow(models.parents),]<- -models.parents[2:nrow(models.parents),]
	ind<-which(score.causal!=Inf)
	score.causal[ind]<-score.causal[ind]*(-1)
	
	return (list(score.matrix=score.causal,score.local=models.parents))
}
