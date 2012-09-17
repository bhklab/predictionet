## function to export a network fit with 'netinf' to Cytoscape, possibly with edge-specific and node-specific statistics
## object: object returns by 'netinf' or 'netinf.cv'
## edge.info: matrix of values representing the statistics for each edge; parents in rows, children in columns. A list of matrices could be provided, names of the list will then be used to describe the statistics in cytoscape
## node.info: vector of values representing the statistics for each node; parents in rows, children in columns. A list of vectors could be provided, names of the list will then be used to describe the statistics in cytoscape
## return a igraph object and create a GML file
'netinf2gml' <- 
function(object, edge.info, node.info, file="predictionet") {
	#require(igraph)
	## adjacency matrix representing the topology; parents in rows, children in columns
	if(!missing(edge.info)) { 	edge.info.new <- edge.info }else{ edge.info.new <- NULL}
	net.topo <- object$topology
	if(any(dim(net.topo) <= 0)) { stop("network should contain at least one node!") }
	## matrix of coeffcients for regrnet
	net.topo.coeff <- object$topology.coeff
	if(all(c("prediction.score.cv", "edge.stability", "edge.relevance") %in% names(object))) {
		## object previously computed by 'netinf.cv'
		edge.info <- list("relevance"=object$edge.relevance, "stability"=object$edge.stability)
		node.info <- lapply(object$prediction.score.cv, function(x) { return(apply(X=x, MARGIN=2, FUN=mean, na.rm=FALSE)) })
	}

	
	if(!is.null(edge.info.new)) {
		if(!is.list(edge.info.new)) { edge.info.new <- list("statistic"=edge.info.new)}

		edge.info.new <- edge.info.new[sapply(edge.info.new, function(x) { return(!is.null(x) && !all(is.na(x))) })]	
		## check dimensions of edge.info
		if(!all(sapply(edge.info.new, dim) == nrow(net.topo))) { stop("edge.info should be a matrix or a list of matrices with the same dimensions that the adjacency describing the network topology!") }else{
		edge.info <- c(edge.info,edge.info.new)	
		}
	}
	if(!missing(node.info)) {
		if(!is.list(node.info)) { node.info <- list("statistic"=node.info) }
		node.info <- node.info[sapply(node.info, function(x) { return(!is.null(x) && !all(is.na(x))) })]
		## check dimensions of node.info
		if(!all(sapply(node.info, length) == nrow(net.topo))) { stop("node.info should be a vector or a list of vectors of length equal to the number of rows/columns of the adjacency describing the network topology!") }
	}

	## create igraph object
	net.igraph <- igraph0::graph.adjacency(adjmatrix=net.topo, mode="directed")
	## graph attributes
	net.igraph <- igraph0::set.graph.attribute(graph=net.igraph, name="Source", value=sprintf("predictionet (R) package version %s", sessionInfo(package="predictionet")$otherPkgs[[1]]$Version))
	## edge attributes
	ee <- igraph0::E(graph=net.igraph)
	## colors for the edges
	rr <- rep("#000000", length(ee))
	net.igraph <- igraph0::set.edge.attribute(graph=net.igraph, name="color", index=ee, value=rr)
	## edge info
	if(!missing(edge.info) && ecount(net.igraph) > 0) {
		for(i in 1:length(edge.info)) {
			if(!is.null(edge.info[[i]])) {
				rr <- edge.info[[i]][igraph0::get.edges(graph=net.igraph, es=ee)+1]
				net.igraph <- igraph0::set.edge.attribute(graph=net.igraph, name=names(edge.info)[[i]], index=ee, value=rr)
			}
		}
		ein <- names(edge.info)
	} else { ein <- NULL }
	## node attributes
	vv <- igraph0::V(graph=net.igraph)
	## color
	rr <- rep("#0099ff", length(vv))
	net.igraph <- igraph0::set.vertex.attribute(graph=net.igraph, name="color", index=vv, value=rr)
	## node info
	if(!missing(node.info) && vcount(net.igraph) > 0) {
		for(i in 1:length(node.info)) {
			if(!is.null(node.info[[i]])) {
				rr <- node.info[[i]][vv+1]
				net.igraph <- igraph0::set.vertex.attribute(graph=net.igraph, name=names(node.info)[[i]], index=vv, value=rr)
			}
		}
		nin <- names(node.info)
	} else { nin <- NULL }
	.exportGML(graph=net.igraph, edge.attributes=ein, vertex.attributes=nin, file=sprintf("%s.gml", file))
	## copy Vizmap Propoerty file for cytoscape
	file.copy(from=system.file(file.path("extdata", "preditionet_vizmap2.props"), package="predictionet"), to=file)
	invisible(net.igraph)
}
