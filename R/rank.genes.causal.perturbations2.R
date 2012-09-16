### function ranking the genes based on mrmr+prior+causality for each gene as a target
## data: Matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
## perturbations: matrix of {0, 1} specifying whether a gene has been pertubed in some experiments. Dimensions should be the same than data
## priors: matrix of prior information avalable for gene-gene interaction (parents in columns, children in rows). Values may be probabilities or any other weights (citations count for instance)
## predn: indices of genes to fit during network inference. If missing, all the genes will be used for network inference
## priors.weight: real value in [0, 1] specifying the weight to put on the priors (0=only the data are used, 1=only the priors are used to infer the topology of the network). If 'priors.weight' is missing it will be optimized gene by hene in an automatic way.
## maxparents: maximum number of parents allowed for each gene
### returns matrix: each row corresponds to a target gene, columns contain the ranked genes (non NA values; gene names)
`.rank.genes.causal.perturbations2` <- 
function(priors, data, perturbations, predn, priors.weight, maxparents, causal=TRUE) {
	
	########################
	### initialize variables
	########################
	if(missing(perturbations) || is.null(perturbations)) {
		perturbations <- matrix(FALSE, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	}
	res <- matrix(NA, nrow=length(predn), ncol=ncol(data) - 1, dimnames=list(dimnames(data)[[2]][predn], seq(1, ncol(data) - 1)))
	res.netw <- matrix(-1, nrow=ncol(data), ncol=length(predn), dimnames=list(dimnames(data)[[2]], dimnames(data)[[2]][predn]))
	data.original <- data.matrix(data)
	score.causal <- matrix(Inf, nrow=ncol(data), ncol=ncol(data), dimnames=list(colnames(data), colnames(data)))
	if(priors.weight != 1) {
		## data will be used to rank the potential parent genes
		mrmr.global <- NULL	
		## no perturbations
		## compute mrmr score for the target genes
		
#		if(causal){
#			local.maxparents<-floor(min(max(3*maxparents,0.05*ncol(data)),0.5*(ncol(data))))
#		}else{
#			local.maxparents<-maxparents
#		}
		
		
		mrmr.global <- .Call("mrnet_adapted2", data.original,as.integer(is.na(data.original)), priors, priors.weight, maxparents, ncol(data.original), nrow(data.original), predn, length(predn), -1000)
		mrmr.global <- matrix(mrmr.global, nrow=ncol(data.original), ncol=ncol(data.original), byrow=FALSE)
		dimnames(mrmr.global) <- list(colnames(data.original), colnames(data.original))

		mrmr.global[mrmr.global == -1000] <- NA


		for(i in 1:length(predn)) {
			data <- data.original
			## remove observation where the target variable was perturbed (no link between observation and target variable)
			ind <- which(is.na(perturbations[, predn[i]]) | perturbations[, predn[i]] == 1)
			if(length(ind)>0) {
				data <- data[-ind,]
				## since the target gene has been perturbed in some of the experiments we should recompute the mim and mrmr matrices
				## compute mrmr score for the target genes
				## mrmr_adapted(SEXP data, SEXP maxparents, SEXP nvar, SEXP nsample, SEXP predn, SEXP npredn, SEXP threshold);
				mrmr <- .Call("mrnet_adapted2", data,as.integer(is.na(data)), priors, priors.weight, maxparents, ncol(data), nrow(data), predn, length(predn), -1000)
				mrmr <- matrix(mrmr, nrow=ncol(data), ncol=ncol(data), byrow=FALSE)
				dimnames(mrmr) <- list(colnames(data),colnames(data))
				mrmr[mrmr == -1000] <- NA
			} else { mrmr <- mrmr.global}
			## mrmr scores on the diagonal should be NA
			diag(mrmr) <- NA

			## select only the variables specified by the user (predn)
			mrmr <- mrmr[ , predn, drop=FALSE]
			########################
			### compute mrmr score matrix
			########################
			## normalization of mrmr scores
			mrmr <- mrmr/max(abs(mrmr), na.rm=TRUE)
			## avoid mrmr score of exactly 0
			mrmr[!is.na(mrmr) & mrmr == 0] <- 1e-16
			score <- matrix(0, nrow=ncol(data), ncol=ncol(data), dimnames=list(colnames(data), colnames(data)))
			score[ , predn][!is.na(mrmr)] <- mrmr[!is.na(mrmr)]
			
			if(causal) {
				########################
				## find all adjacencies of node to predict in score
				########################
				tmp.adj <- matrix(0,ncol=ncol(data),nrow=ncol(data), dimnames=list(colnames(data), colnames(data)))
				tmp.adj[predn[i],] <- tmp.adj[ ,predn[i]] <- score[ ,predn[i]]
				
				########################
				## determine all triplets X.k-X.j-X.l from these adjecencies
				########################
				trip <- (.network2triplets(tmp.adj,predn[i]))
				parents <- NULL
				score.parents <- NULL
				if(length(trip)>0){

					########################
					## compute I(X.k,X.l)-I(X.k,X.l|X.j)
					########################
					trip <- as.data.frame(cbind(trip,(.get.ii4triplets.gaussian(data,trip))))
					########################
					## remove those triplets for which the outer nodes are connected
					########################
					ind.rm <- NULL
					for(j in 1:nrow(trip)){
						if(score[trip[j,1],trip[j,3]] !=0 || score[trip[j,3],trip[j,1]] !=0) {
							ind.rm <- c(ind.rm,j)
						}
					}
					if(!is.null(ind.rm)) { trip <- (trip[-ind.rm, , drop=FALSE]) }
					########################
					## if there are triplets for which I(X,Z)-I(X,Z|Y)<0, add these to the list of parents
					########################
					if(nrow(trip)>0){
						for(j in 1:nrow(trip)){
							if(!is.na(trip[j,4]) && !is.nan(trip[j,4]) && trip[j,4] < score.causal[trip[j,1],trip[j,2]]){
								score.causal[trip[j,1],trip[j,2]] <- trip[j,4]
							}
							if(!is.na(trip[j,4]) && !is.nan(trip[j,4]) && trip[j,4] < score.causal[trip[j,3],trip[j,2]]){
								score.causal[trip[j,3],trip[j,2]] <- trip[j,4]
							}
						}
					}
				}
			} else {
			## no causality inference, all connected nodes are parents
				res.netw[ ,i] <- mrmr[ , predn[i]]
				res.netw[is.na(res.netw[ ,i]),i] <- -1
			}
		}
	}

	if(causal) {
		if(!all(is.infinite(score.causal))) { score.causal <- score.causal/(-max(abs(score.causal[!is.infinite(score.causal)]))) }
		score.causal[is.na(score.causal) | is.infinite(score.causal)] <- -1
		res.netw <- score.causal[,predn]
	}
	

	########################
	## (1-w)*M+w*P weighted score from prior knowledge and data
	########################
	score <- (1-priors.weight)*res.netw+priors.weight*priors[ , predn, drop=FALSE]
	## no self-loops
	diag(score[dimnames(score)[[2]], dimnames(score)[[2]]]) <- -1
	for(j in 1:length(predn)){
		tmp.s <- sort(score[,j],decreasing=TRUE,index.return=TRUE)
		ind.rm <- (which(tmp.s$x<=0))
		if(length(tmp.s$x[-ind.rm])>0){
			res[j,(1:(length(tmp.s$x[-ind.rm])))] <- names(tmp.s$x[-ind.rm])
			if(length(tmp.s$x[-ind.rm])>maxparents){
				res[j,((maxparents+1):(length(tmp.s$x[-ind.rm])))] <- NA
			}
		}
	}
	return(list("parents"=res, "score"=score))
}
