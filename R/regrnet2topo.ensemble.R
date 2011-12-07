## Function transforming an object received from `.fit.regrnet.causal.ensemble` into an adjacency matrix
## containing either 1 for each inferred edge or the edge.stability value for that edge 
## Note that the ensemble  topology has target variables as columns thus allowing for multiple columns the same name
## net: an object returned from `.fit.regrnet.causal.ensemble`
## coefficients: default=FALSE returning an adjacency matrix with 1 for each inferred edge, =TRUE returns the corresponding edge.relevance for each inferred edge

`.regrnet2topo.ensemble`<- function(net, coefficients=FALSE) {
	geneid <- (net$input)
	nr <- length(geneid)
	
	res <- matrix(0, nrow=length(net$varnames), ncol=length(net$input), dimnames=list(net$varnames,net$input))
	for(i in 1:length(net$input)){
		res[,i]<-net$edge.relevance[[i]]
	}

	if(!coefficients){
		res[res!=0]<-1
	}
	return(res)
}
