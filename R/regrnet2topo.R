### function transforming a regression model into an adjacency matrix
## net: regression model for each gene
## coeffcients: if TRUE, the coefficients are extracted from the regression models for each existing edges,  boolean representing the presence of the edges otherwise
### returns adjacency matrix: parents in rows and children in columns (res[i,j]==1 means edge from i to j)
`.regrnet2topo` <- 
function(net, coefficients=FALSE) {
	geneid <- names(net$input)
	nr <- length(geneid)
	net.model <- net$model
	res <- matrix(0, nrow=length(net$varnames), ncol=nr, dimnames=list(net$varnames,geneid))
	beta_0<-NULL
	for (i in 1:nr) {
		model.i <- net.model[[i]]
		if(!is.numeric(model.i)) { ## non null model
			if(coefficients) {
				beta <- coefficients(object=model.i)
				beta_0 <- c(beta_0,beta[1])
				res[names(beta)[-1], geneid[i]] <- beta[-1]
			} else { res[net$input[[i]], geneid[i]] <- 1 }
		}else{
			if(coefficients) {
				beta_0 <- c(beta_0,model.i)
			}
		}
	}
	
	if(coefficients){
		names(beta_0)<-colnames(res)
		res<-rbind(beta_0,res)
	}
	return(res)
}
