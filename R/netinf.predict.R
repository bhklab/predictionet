### Function predicting values of nodes from a previously inferred network
`netinf.predict` <- 
function(net, data, categories, perturbations, subset, predn) {
	res <- matrix(NA, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	switch(net$method, 
		   "bayesnet"={
		   res <- .pred.onegene.bayesnet.fs(net=net$net, data=data, categories=categories, perturbations=perturbations, subset=subset, predn=predn)
		   }, 
		   "regrnet"={
		   res <- .pred.onegene.regrnet.fs(net=net$net, data=data, perturbations=perturbations, subset=subset, predn=predn)
		   }, 
		   "regrnet.ensemble"={
			## nsemble network model
		   stop("Ensemble regression-based network inference is not implemented yet!")
		   })
	return(res)
}