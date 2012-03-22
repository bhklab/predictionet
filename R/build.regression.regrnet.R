### Function fitting a regression model for each gene in the data
## data: Matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## predn: indices or names of genes to fit during network inference. If missing, all the genes will be used for network inference
## topo: inferred topology (matrix containing {0,1})
### returns a regression network 
`.build.regression.regrnet` <- 
function(net, data, predn, perturbations, regrmodel=c("linear", "linear.penalized"), seed) {
	if(!missing(seed)) { set.seed(seed) }
	topo <- net$topology
	ensemble <- net$ensemble
	if(ensemble) {
		nn <- colnames(topo) ## variables (genes) to fit during network inference
	} else {
		if(missing(predn)) { predn <- seq(1,ncol(data)) } else if(length(intersect(predn, colnames(data))) > 0) { predn <- match(predn,colnames(data)) }
		nn <- dimnames(data)[[2]][predn] ## variables (genes) to fit during network inference
	}
	regrmodel <- match.arg(regrmodel)
	regrnet <- NULL
	########################
	### data preparation
	########################
	dd <- data.frame(data)
		
	########################
	### build regression model for each gene
	########################
	for(i in 1:length(nn)) {
	### consider only those genes for the model which are ranked (not the NA entries)
#	modelv <- myparents[[nn[i]]]
		if(ensemble) {
			modelv <- names(which(topo[,i]>0))
		} else {
			modelv <- names(which(topo[,nn[[i]]]>0))
		}
		if(length(modelv)==0) { ## none variables have been selected
			mm.reg <- mean(dd[ , nn[i]], na.rm=TRUE) ## null model (only intercept)
		} else {
			modelv <- sort(modelv) ## sort variables by name and remove NA
			ff <- formula(sprintf("%s ~ 1 + %s", nn[i], paste(modelv, collapse=" + ")))
			dd2 <- dd
			dd2[perturbations[ , nn[i]], nn[i]] <- NA
			switch(regrmodel, 
				   "linear"={
				   mm.reg <- lm(formula=ff, data=dd2)
				   },
				   "linear.penalized"={
				   #require(penalized)
					ccix <- !perturbations[ , nn[i]]
				   optlambda1 <- penalized::optL1(response=ff, data=dd2[ccix, ,drop=FALSE], model="linear", lambda2=0, minlambda1=1, maxlambda1=10, trace=FALSE, fold=10)
				   mm.reg <- penalized::penalized(response=ff, data=dd2[ccix, ,drop=FALSE], model="linear", lambda1=optlambda1$lambda, lambda2=0, trace=FALSE)
				   })
		}
		regrnet <- c(regrnet, list(mm.reg))
	}
	names(regrnet) <- nn
	return(c(net, "lrm"=list(.regrnet2matrixtopo(c(net, "lrm"=list(regrnet))))))

#	return(c(net, "lrm"=list(regrnet), "topo.coeff"=list(.regrnet2matrixtopo(c(net, "lrm"=list(regrnet))))))
}
