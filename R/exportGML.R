## inspire from Gábor Csárdi
## source: http://tolstoy.newcastle.edu.au/R/e13/help/11/04/10607.html
`.exportGML` <- 
function(graph, edge.attributes, vertex.attributes, file) {
	#require(igraph)
	if(!missing(edge.attributes) && !is.null(edge.attributes) && !all(edge.attributes %in% igraph::list.edge.attributes(graph))) { stop("edge attributes are not stored in the igraph object!") }
	if(!missing(vertex.attributes) && !is.null(vertex.attributes) && !all(vertex.attributes %in% igraph::list.vertex.attributes(graph))) { stop("vertex attributes are not stored in the igraph object!") }
	
	file <- file(file, "w")
	## graph attributes
	cat("Creator \"igraph exportGML\"\n", file=file)
	cat("Version 1.0\n", file=file)
	cat("graph\t[\n", file=file)
	cat("\tdirected", as.integer(igraph::is.directed(graph)), "\n", file=file)
	
	## vertex attributes
	if(vcount(graph) > 0) {
		for (i in seq_len(igraph::vcount(graph))) {
			cat("\tnode\t[\n", file=file)
			cat("\t\tid", i-1, "\n", file=file)
			cat("\t\tgraphics\t[\n", file=file)
			cat("\t\t\tfill\t\"", igraph::V(graph)$color[i], "\"\n", sep="", file=file)
			cat("\t\t\ttype\t\"ellipse\"\n", file=file)
			cat("\t\t\toutline\t\"#cccccc\"\n", file=file)
			cat("\t\t\toutline_width\t1.0\n", file=file)
			cat("\t\t]\n", file=file)
			cat("\t\tlabel\t\"", igraph::V(graph)$name[i], "\"\n", file=file)
			## special vertex attributes
			if(!missing(vertex.attributes) && !is.null(vertex.attributes)) {
				for(jj in 1:length(vertex.attributes)) {
					aa <- igraph::get.vertex.attribute(graph, name=vertex.attributes[jj])
					if(is.numeric(aa)) {
						#cat("\t\t", edge.attributes[jj], "\t", sprintf("%.3f", aa[i]), "\n", file=file)
						cat(sprintf("\t\t%s\t%.3f\n", vertex.attributes[jj], aa[i]), file=file)
					} else { cat(sprintf("\t\t%s\t\"%s\"\n", vertex.attributes[jj], aa[i]), file=file) }
				}
			}
			cat("\t]\n", file=file)
		}
	}
	
	## edge attributes
	if(ecount(graph) > 0) {
		el <- igraph::get.edgelist(graph, names=FALSE)
		eln <- apply(igraph::get.edgelist(graph, names=TRUE), 1, paste, collapse="-")
		for (i in seq_len(nrow(el))) {
			cat("\tedge\t[\n", file=file)
			#cat("\t\tlabel\t", eln[i], "\n", file=file)
			cat("\t\tsource", el[i,1], "\n", file=file)
			cat("\t\ttarget", el[i,2], "\n", file=file)
			cat("\t\tgraphics\t[\n", file=file)
			cat("\t\t\twidth\t1.0\n", file=file)
			cat("\t\t\tfill\t\"", igraph::E(graph)$color[i], "\"\n", sep="", file=file)
			cat("\t\t\ttype\t\"line\"\n", file=file)
			cat("\t\t\tsource_arrow\t0\n", file=file)
			cat("\t\t\ttarget_arrow\t3\n", file=file)
			cat("\t\t]\n", file=file)
			## special edge attributes
			if(!missing(edge.attributes) && !is.null(edge.attributes)) {
				for(jj in 1:length(edge.attributes)) {
					aa <- igraph::get.edge.attribute(graph, name=edge.attributes[jj])
					if(is.numeric(aa)) {
						#cat("\t\t", edge.attributes[jj], "\t", sprintf("%.3f", aa[i]), "\n", file=file)
						cat(sprintf("\t\t%s\t%.3f\n", edge.attributes[jj], aa[i]), file=file)
					} else { cat(sprintf("\t\t%s\t\"%s\"\n", edge.attributes[jj], aa[i]), file=file) }
				}
			}
			cat("\t]\n", file=file)
		}
	}
	cat("]\n", file=file) 
	close(file) 
}
