`.list2matrixens` <- 
function(myres) {
	res.edge.relevance<-matrix(0,nc=ncol(myres$topology),nr=nrow(myres$topology),dimnames=dimnames(myres$topology))


	for(i in 1:ncol(myres$topology)){
		res.edge.relevance[,i]<-myres$edge.relevance[[i]]
	}
	

	return(list("method"=myres$method,"topology"=myres$topology,"edge.relevance"=res.edge.relevance))
}
