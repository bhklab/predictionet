### function transforming a regression model into an adjacency matrix
## net: bayesian network
### returns adjacency matrix: parents in rows and children in columns (res[i,j]==1 means edge from i to j)
`.bayesnet2topo` <- 
function(net, coefficients=FALSE) {
	if(coefficients) { return(NULL) }
	#require(catnet)
	return(t(catnet::cnMatParents(net$model)))
}
