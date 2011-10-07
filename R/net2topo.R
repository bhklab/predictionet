### function transforming a regression model into an adjacency matrix
## net: regression model for each gene
## coeffcients: if TRUE, the coefficients are extracted from the regression models for each existing edges,  boolean representing the presence of the edges otherwise
### returns adjacency matrix: parents in rows and children in columns (res[i,j]==1 means edge from i to j)
`net2topo` <- 
function(net, coefficients=FALSE) {
	switch(net$method, 
		"regrnet"={
		   res <- .regrnet2topo(net=net$net, coefficients=coefficients)
		},
		"bayesnet"={
		   res <- .bayesnet2topo(net=net$net, coefficients=coefficients)
		}, 
		{
		   stop("method used to infer the network is unknown!")
		})
	return(res)
}
