`.list2matrixens` <- 
function(netmp) {
	res.edge.relevance <- matrix(0, ncol=ncol(netmp$topology), nrow=nrow(netmp$topology), dimnames=dimnames(netmp$topology))
	for(i in 1:ncol(netmp$topology)){
		res.edge.relevance[,i] <- netmp$edge.relevance[[i]]
	}
	netmp$edge.relevance <- res.edge.relevance
	return(netmp)
}
