## function to identify all children of a parent
## adjmat: adjacency matrix; parents in rows, children in columns
## from: name or index of the parent node to start from
## depth: how deep eth children should be identified
## return a two-column matrix containing the names of the children in the first column and their corresponding depth in the descent in the second column
'adj.descent' <- 
function(adjmat, from, depth) {
	
	####################
	## internal functions
	####################
	## adjmat: adjacency matrix
	## nodest: status of nodes (vector); 0 if not visited, 1 otherwise
	## nodix.from: index of the previous node
	## nodix: index of the current node
	'.adj.DFS' <- 
	function(adjmat, nodix.from, nodix, nodest, kids, depth=0, depth.max) {
		if(nodix.from %in% (1:length(nodest)) && nodest[nodix] == 1) { ## the current node has already been visited
			return(list("kids"=c(NA, depth), "nodest"=nodest))
		} else {
			## the current node has never been visited
			nodest[nodix] <- 1
			kids <- c(nodix, depth)
			ee <- which(adjmat[nodix, ] > 0) ## all children of nodix
			if(length(ee) > 0 && depth < depth.max) { ## some children to look at
				for(i in 1:length(ee)) { ## for all children
					rr <- .adj.DFS(adjmat=adjmat, nodix.from=nodix, nodix=ee[i], nodest=nodest, depth=depth+1, depth.max=depth.max)
					kids <- rbind(kids, rr$kids)
					nodest <- rr$nodest
				}
			}
			return(list("kids"=kids, "nodest"=nodest))
		}
	}
	####################
	if(!is.matrix(adjmat) && nrow(adjmat) != ncol(adjmat)) { stop("the adjacency matrix should be square!") }
	if(is.null(rownames(adjmat))) { rownames(adjmat) <- colnames(adjmat) <- paste("X", 1:nrow(adjmat), sep="") }
	nNames <- rownames(adjmat)
	adjmat <- adjmat[nNames, nNames, drop=FALSE]
	if(missing(depth) || (depth <= 0 || depth >= ncol(adjmat))) { depth <- nrow(adjmat) - 1 }
	if(is.character(from)) { from <- match(from, nNames)}
	
	nodest <- rep(0, ncol(adjmat))
	names(nodest) <- colnames(adjmat)
	#nodest[from] <- 1
	rr <- .adj.DFS(adjmat=adjmat, nodix.from=-1, nodix=from, nodest=nodest, kids=NULL, depth=0, depth.max=depth)
	kids <- rr$kids
	if(!is.matrix(kids)) { return(NULL) } else {
		kids <- kids[-1, , drop=FALSE]
		dimnames(kids) <- list(NULL, c("children", "depth"))
		kids <- kids[order(kids[ , 2], decreasing=FALSE), , drop=FALSE]
		kids <- kids[!is.na(kids[ , 1]) & !duplicated(kids[ , 1]), , drop=FALSE]
		kids[ , 1] <- nNames[kids[ , 1]]
		return(kids)
	}
}
