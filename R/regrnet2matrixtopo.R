### function transforming a regression model into an adjacency matrix
## topo: inferred topology, adjacency matrix containing {0,1}
## net.model: regression model for each gene
### returns adjacency matrix with entries corresponding to the regression coefficients: parents in rows and children in columns (res[i,j]==1 means edge from i to j)
`.regrnet2matrixtopo` <- 
function(net) {
	nr <- ncol(net$topology)
	res <- net$topology
	net.model <- net$lrm
	beta_0 <- NULL
	for (i in 1:nr) {
		model.i <- net.model[[i]]
		if(!is.numeric(model.i)) { ## non null model
				beta <- coefficients(object=model.i)
				beta_0 <- c(beta_0,beta[1])
				res[names(beta)[-1], i ] <- beta[-1]
		}else{
				beta_0 <- c(beta_0,model.i)
		}
	}
	names(beta_0) <- colnames(res)
	res <- rbind(beta_0,res)
	return(res)
}
