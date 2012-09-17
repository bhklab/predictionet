### Function computing the press statistic for all target variables in topology
## topo: inferred topology
## data: matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## returns press statistic for all target variables

`predictionet.press.statistic` <- function(topo, data, ensemble=FALSE, perturbations=NULL) {

	if(missing(perturbations) || is.null(perturbations)) {
		perturbations <- matrix(FALSE, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	} else {
		if(nrow(perturbations) == 1) {
			perturbations[1, ] <- as.logical(perturbations[1, ])
		} else { perturbations <- apply(perturbations, 2, as.logical) }
		dimnames(perturbations) <- dimnames(data)
	}
	
	if(ensemble){
		mypert.ens <- NULL
		for(i in 1:ncol(topo)){
			mypert.ens <- cbind(mypert.ens,perturbations[,colnames(topo)[i]])
		}
		colnames(mypert.ens) <- colnames(topo)
		perturbations <- mypert.ens
	}
	
	res <- matrix(0,ncol=ncol(topo),nrow=nrow(data),dimnames=list(rownames(data),colnames(topo)))
	vec.nsamples <- rep(0,ncol(topo))
	
	for(i in 1:ncol(topo)){
			target <- colnames(topo)[i]
			ind <- which(topo[,i]!=0)
			ind.pert <- which(perturbations[,i]==1)
			vec.nsamples[i] <- nrow(data)-length(ind.pert)
			if(length(ind.pert)>0) {
				mydata <- data.matrix(data[-ind.pert, ])
			} else {
				mydata <- data.matrix(data)
			}
			if(length(ind)>0) {
				if(length(ind)==1) {
					res[rownames(mydata),i] <- .regrlin(mydata[,ind, drop=FALSE], mydata[,target])
				} else {
					res[rownames(mydata),i] <- .regrlin(mydata[,ind, drop=FALSE], mydata[,target])
				}
			} else {
				res[rownames(mydata),i] <- .regrlin(rep(1,nrow(mydata)), mydata[,target])
			}
	}
	res <- res^2
	res <- apply(res, 2, sum)
	res <- res/vec.nsamples
	return(res)
}
