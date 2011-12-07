### Function predicting values of nodes from a previously inferred network
## net: object returned by netinf 
## data: matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## categories: Either an integer specifying the number of categories used to discretize the variables or a list of categories for each of the variables in 'data' matrix. Categories should be integers from 1 to n. If method='bayesnet', this parameter should be specified by the user.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## subset: vector of indices to select only subset of the observations
## predn: indices or names of genes to fit during network inference. If missing, all the genes will be used for network inference
## returns matrix containg a prediction for each sample and each variable in predn

`netinf.predict` <- 
function(net, data, categories, perturbations, subset, predn, ensemble=FALSE, topo.coeff=FALSE) {
	res <- matrix(NA, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	switch(net$method, 
		   "bayesnet"={
			res <- .pred.onegene.bayesnet.fs(net=net$net, data=data, categories=categories, perturbations=perturbations, subset=subset, predn=predn)
		   }, 
		   "regrnet"={
		   if(ensemble){
				if(missing(predn)){predn<-colnames(net$topology)}else if(length(intersect(predn,colnames(net$topology)))>0){predn<-match(predn,colnames(data))}
		   }else{
				if(missing(predn)){predn<-seq(1,ncol(data))}else if(length(intersect(predn,colnames(data)))>0){predn<-match(predn,colnames(data))}
		   }
		   if(topo.coeff==FALSE && ensemble){
				stop("The input options don't agree on what to do, set topo.coeff=TRUE if you have an ensemble object, this needs to contain the topology.coeff information")
		   }
		   if(topo.coeff && ensemble){
				res <- .pred.onegene.regrnet.fs(topo.coeff=net$topology.coeff, data=data, perturbations=perturbations, subset=subset, predn=predn, ensemble=ensemble)
		   }else{
				res <- .pred.onegene.regrnet.fs(topo.coeff=.regrnet2matrixtopo(net$topology[,predn,drop=FALSE] , build.regression.regrnet(data,predn,net$topology)), data=data, perturbations=perturbations, subset=subset, predn=predn, ensemble=ensemble)
		   }
		   }
	)
	return(res)
}
