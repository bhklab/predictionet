### function transforming a regression model into an adjacency matrix
## net: regression model for each gene
## regr.model: if adjacency matrix should be filled with regression coefficients, the regression model has to be provided
## coeffcients: if TRUE, the coefficients are extracted from the regression models for each existing edges,  boolean representing the presence of the edges otherwise
### returns adjacency matrix: parents in rows and children in columns (res[i,j]==1 means edge from i to j)
`net2topo` <- 
function(net, regr.model=NULL, coefficients=FALSE) {
	switch(net$method, 
		"regrnet"={
		   if(coefficients){
		   res <- .regrnet2matrixtopo(topo=net$topology, regr.model)
		   }else{
		   res<-net$topology
		   }
		},
		"bayesnet"={
		   res <- .bayesnet2topo(net=net$net, coefficients=coefficients)
		}, 
		{
		   stop("method used to infer the network is unknown!")
		})
	return(res)
}
