`.topo2stab`<-function(topo){
	tmp.topo<-topo
	tmp.names<-intersect(colnames(topo[[1]]),colnames(topo[[1]]))
	res<-matrix(0,nc=length(tmp.names),nr=nrow(topo[[1]]),dimnames=list(rownames(topo[[1]]),tmp.names))
	counts<-rep(0,length(tmp.names))
	names(counts)<-tmp.names
	
	for(j in 1:length(tmp.topo)){
		topo<-tmp.topo[[j]]
		
		for(i in 1:ncol(topo)){
			counts[colnames(topo)[i]]<-counts[colnames(topo)[i]]+1
			res[which(topo[,i]==1),colnames(topo)[i]]<-res[which(topo[,i]==1),colnames(topo)[i]]+1
		}
	}
	for(i in 1:ncol(res)){
		res[,i]<-res[,i]/counts[i]
	}
	return(res)
}
