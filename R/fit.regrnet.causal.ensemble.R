`.fit.regrnet.causal.ensemble` <- 
function(res.causal,models, data, perturbations, priors, predn, maxparents=3, priors.weight=0.5, causal=TRUE, seed) {
	if(causal && maxparents > (ncol(data) * 0.5)) { warning("maximum number of parents may be too large, causal inference requires sparsity in the inferred network; please decrease maxparents parameter for better results!") }
	if(!missing(seed)) { set.seed(seed) }
	if(missing(predn) || is.null(predn) || length(predn) == 0) { predn <- 1:ncol(data) } else { if(is.character(predn)) { predn <- match(predn, dimnames(data)[[2]]) } else { if(!is.numeric(predn) || !all(predn %in% 1:ncol(data))) { stop("parameter 'predn' should contain either the names or the indices of the variables to predict!")} } }
	
	regrnet <- NULL
########################
### data preparation
########################
	dd <- data.frame(data)
	diag(priors) <- 0 ## no self-loops
	
#	res.causal <- rank.genes.causal.ensemble(res.ensemble,dd)
#	models <- res.ensemble$models
	res.causal.local <- res.causal$score.local
	nn <- models[1,]
	
	edge.relevance <- NULL

	for(i in 1:ncol(res.causal.local)){
		tmp.causal <- (1-priors.weight)*res.causal.local[,i]
		names(tmp.causal) <- models[,i]
		if(priors.weight>0 && priors.weight<1){
			ind.prior <- which(priors[,models[1,i]]!=0)
			if(length(ind.prior)>0){
				for(j in 1:length(ind.prior)){
					ind.replace <- which(names(tmp.causal)==names(ind.prior[j]))
					if(length(ind.replace)>0){
						tmp.causal[ind.replace] <- tmp.causal[ind.replace]+priors.weight*priors[ind.prior[j],models[1,i]]
					}else{
						tmp.names <- names(tmp.causal)
						tmp.causal <- c(tmp.causal,priors.weight*priors[ind.prior[j],models[1,i]])
						names(tmp.causal) <- c(tmp.names,names(ind.prior[j]))
					}
				}
			}
		}else if(priors.weight==1){
			ind.prior <- which(priors[,models[1,i]]!=0)
			if(length(ind.prior)>0){
				tmp.causal <- c(models[1,i],priors[ind.prior,models[1,i]])
				names(tmp.causal) <- c(models[1,i],names(ind.prior))
			}else{
				tmp.causal <- NULL
			}
		}
		
		tmp.causal <- tmp.causal[tmp.causal>0]
		tmp.relevance <- rep(0,ncol(data))
		names(tmp.relevance) <- colnames(data)
		
		if(!is.null(tmp.causal) && length(tmp.causal)>1){
			tmp.causal <- tmp.causal[-1]
			tmp.causal <- (sort(tmp.causal,decreasing=TRUE))
			
#			if(length(tmp.causal)>0){
#		modelv <- sort(names(tmp.causal[1:min(maxparents,length(tmp.causal))]))
#				ff <- formula(sprintf("%s ~ 1 + %s", nn[i], paste(modelv, collapse=" + ")))
				
				tmp.relevance[names(tmp.causal[1:min(maxparents,length(tmp.causal))])] <- tmp.causal[1:min(maxparents,length(tmp.causal))]
#			}else{
#				ff <- formula(sprintf("%s ~ 1", nn[i]))
#			}
		}
#else{
#			ff <- formula(sprintf("%s ~ 1", nn[i]))
#		}
#		mm.reg <- lm(formula=ff, data=dd)
#		regrnet <- c(regrnet, list(mm.reg))
		edge.relevance <- c(edge.relevance,list(tmp.relevance))

	}
#	names(regrnet) <- nn
	names(edge.relevance) <- nn
#return(list("input"=models[1,],"varnames"=colnames(data), "model"=regrnet, "edge.relevance"=edge.relevance))
return(list("input"=models[1,],"varnames"=colnames(data), "edge.relevance"=edge.relevance))
	
}
