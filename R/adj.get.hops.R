## function to compute the shortest path
'adj.get.hops' <- 
function(adjmat) {
	#require(RBGL)
	hops <- RBGL::floyd.warshall.all.pairs.sp(as(adjmat,"graphNEL"))
	#dimnames(hops) <- dimnames(adjmat)
	hops[abs(hops) == Inf] <- NA
	return(hops)
}
