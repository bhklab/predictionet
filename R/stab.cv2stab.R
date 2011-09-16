`.stab.cv2stab`<-function(stab.cv, topo){
		tmp.names<-intersect(colnames(topo),colnames(topo))
		
		res<-matrix(0,nc=length(tmp.names),nr=nrow(topo),dimnames=list(rownames(topo),tmp.names))
		stab<-res
		for(i in 1:ncol(topo)){
			res[which(topo[,i]==1),colnames(topo)[i]]<-res[which(topo[,i]==1),colnames(topo)[i]]+1
		}
		stab[res != 0] <- stab.cv[res != 0]
		return(stab)
	}
