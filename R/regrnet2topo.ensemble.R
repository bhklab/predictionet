`.regrnet2topo.ensemble` <- 
function(net, coefficients=FALSE) {
	geneid <- (net$input)
	nr <- length(geneid)
	res.ensemble <-NULL
	
	for(i in 1:length(net$model)){
		model.i <- net$model[[i]]
		res <- matrix(0, nrow=length(net$varnames), ncol=length(net$varnames), dimnames=list(net$varnames,net$varnames))
		if(!is.numeric(model.i)) { ## non null model
			if(coefficients) {
				beta <- coefficients(object=model.i)
				res[names(beta)[-1], geneid[i]] <- beta[-1]
				res<-rbind(rep(beta[1],ncol(res)),res)
				dimnames(res)<-list(c("beta_0",net$varnames),net$varnames)
			} else { 
				beta <- names(coefficients(object=model.i))
				res[beta[-1], geneid[i]] <- 1 
			}
		}
		
		res.ensemble<-c(res.ensemble,list(res[,geneid[i]]))
	}
	names(res.ensemble)<-geneid
	return(res.ensemble)
}
