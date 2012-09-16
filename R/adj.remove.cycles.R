## function to remove cycles in an adjacency matrix
## adjmat: adjacency matrix; parents in rows, children in collumns
## return the adjacency matrix with the lowest entries being removed such that the network is now acyclic
'adj.remove.cycles' <- 
function(adjmat, from, maxlength) {
pathlength=0	
	####################
	## internal functions
	####################
	## adjmat: adjacency matrix
	## nodest: status of nodes (vector); 0 if not visited, 1 otherwise
	## nodix.from: index of the previous node
	## nodix: index of the current node
	'.adj.remove.cycles.DFS' <- 
	function(adjmat, adjmat.mask, nodix.from, nodix, nodest, pathlength, maxlength) {
		
		if(nodix.from %in% (1:length(nodest)) && nodest[nodix] == 1) { ## the current node has already been visited
			## remove edge
			adjmat.mask[nodix.from, nodix] <- TRUE
			return(adjmat.mask)
		} else {
			if(pathlength > maxlength) {
				return(adjmat.mask)
			} else {
				## the current node has never been visited
				## remove the edges previously identified as responsible for a cycle
				adjmat2 <- adjmat
				adjmat2[adjmat.mask] <- 0
				ee <- which(adjmat2[nodix, ] > 0) ## all children of nodix
				ee <- ee[order(adjmat2[nodix, ee], decreasing=FALSE)] ## order them based on number of citations
				if(length(ee) > 0) { ## some children to look at
					for(i in 1:length(ee)) { ## for all children
						nodest2 <- nodest
						nodest2[nodix] <- 1
						adjmat.mask <- .adj.remove.cycles.DFS(adjmat=adjmat, adjmat.mask=adjmat.mask, nodix.from=nodix, nodix=ee[i], nodest=nodest2, pathlength=pathlength+1,maxlength=maxlength)
					}
					return(adjmat.mask)
				} else { return(adjmat.mask) }
			}
		}
	}
	####################
	if(!is.matrix(adjmat) && nrow(adjmat) != ncol(adjmat)) { stop("the adjacency matrix should be square!") }
	if(is.null(rownames(adjmat))) { rownames(adjmat) <- colnames(adjmat) <- paste("X", 1:nrow(adjmat), sep="") }
	nNames <- rownames(adjmat)
	adjmat <- adjmat[nNames, nNames, drop=FALSE]
	adjmat2 <- adjmat
	adjmat.mask <- matrix(FALSE, nrow=nrow(adjmat), ncol=ncol(adjmat), dimnames=dimnames(adjmat))
	
	## remove bidirectional edges wrt the number of interactions
	## on the diagonal
	iix <- which(diag(adjmat2) != 0)
	if(length(iix) > 0) {
		## some entries in the diagonal should be removed
		adjmat.mask[cbind(iix, iix)] <- TRUE
	}
	adjmat2[adjmat.mask] <- 0
	## bidirectional edges
	#iix <- which(upper.tri(adjmat2), arr.ind=TRUE)
	#iix2 <- t(apply(iix, 1, function(x, y) { if(y[x[1], x[2]] <= y[x[2], x[1]]) { return(x) } else { return(rev(x)) } }, y=adjmat2))
	#iix2 <- iix2[which(adjmat2[iix2] != 0), ]
	#adjmat.mask[iix2] <- TRUE
	#adjmat2[adjmat.mask] <- 0
	
	## order nodes with the maximum prior
	if(missing(from)) {
	## for all the nodes cycles will be removed
		nnix <- order(apply(adjmat2, 1, max), decreasing=FALSE)
	} else {
		if(is.integer(from)) { nnix <- from } else { nnix <- match(from, nNames) }
	}
	for(i in 1:length(nnix)) {
		nodest <- rep(0, ncol(adjmat2))
		names(nodest) <- colnames(adjmat2)
		adjmat.mask <- .adj.remove.cycles.DFS(adjmat=adjmat2, adjmat.mask=adjmat.mask, nodix.from=-1, nodix=nnix[i], nodest=nodest,pathlength=pathlength,maxlength=maxlength)
		adjmat2[adjmat.mask] <- 0
	}
	return(list("adjmat.acyclic"=adjmat2, "adjmat.removed"=adjmat.mask))
}
