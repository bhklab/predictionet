`.list2matrixens` <- 
function(myres) {
	names.target<-names(myres$topology)
	res<-matrix(0,nc=length(myres$topology),nr=length(myres$topology[[1]]),dimnames=list(names(myres$topology[[1]]),names.target))
	res.coeff<-matrix(0,nc=length(myres$topology),nr=length(myres$topology.coeff[[1]]),dimnames=list(names(myres$topology.coeff[[1]]),names.target))
	res.edge.relevance<-res

	for(i in 1:length(myres$topology)){
		res[,i]<-myres$topology[[i]]
		res.coeff[,i]<-myres$topology.coeff[[i]]
		res.edge.relevance[,i]<-myres$edge.relevance[[i]]
	}
	

	return(list("method"=myres$method,"topology"=res,"topology.coeff"=res.coeff,"edge.relevance"=res.edge.relevance))
}
