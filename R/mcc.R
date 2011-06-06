## function to compute the Matthews Correlation Coefficient (MCC) in a classification framework
## ct: confusion table (matrix or table)
## nbcat: the actual number of categories for the truth and the predictions
## Example:
## rr <- sapply(1:1000, function(x, cc) { set.seed(x); ctt <- table(sample(1:cc, 100, replace=TRUE), sample(1:cc, 100, replace=TRUE)); return(mcc(ct=ctt, nbcat=cc)); }, cc=3)
## hist(rr)
`mcc` <- 
function(ct, nbcat=nrow(ct)) {
	if(nrow(ct) != ncol(ct)) { stop("the confusion table should be squared!") }
	if(!(sum(ct)==sum(diag(ct))) &&  (length(which(apply(ct, 1, sum) == 0)) == (nbcat-1) & ((length(which(apply(ct, 2, sum) == 0)) != (nbcat-1)) | (length(which(apply(ct, 2, sum) == 0)) == (nbcat-1)))) || (length(which(apply(ct, 2, sum) == 0)) == (nbcat-1) & ((length(which(apply(ct, 1, sum) == 0)) != (nbcat-1)) | (length(which(apply(ct, 1, sum) == 0)) == (nbcat-1)) & sum(diag(ct)) == 0))) { ct <- ct + matrix(1, nc=nbcat, nr=nbcat) } ### add element to categories if nbcat-1 predictive categories do not contain elements. Not in case where all are correct!
	myx <- matrix(TRUE, nrow=nrow(ct), ncol=ncol(ct))
	diag(myx) <- FALSE
	if(sum(ct[myx]) == 0) { return(1) }
	myperf <- 0
	for(k in 1:nbcat) {
		for(m in 1:nbcat) {
			for(l in 1:nbcat) {
				myperf <- myperf + ((ct[k, k] * ct[m, l]) - (ct[l, k] * ct[k, m]))
			}
		}
	}
	aa <- 0
	for(k in 1:nbcat) {
		cc <- 0
		for(l in 1:nbcat) { cc <- cc + ct[l, k] }
		dd <- 0
		for(f in 1:nbcat) {
			for(g in 1:nbcat) { if(f != k) { dd <- dd + ct[g, f] } }
		}
		aa <- aa + (cc * dd)
	}
	bb <- 0
	for(k in 1:nbcat) {
		cc <- 0
		for(l in 1:nbcat) { cc <- cc + ct[k, l] }
		dd <- 0
		for(f in 1:nbcat) {
			for(g in 1:nbcat) { if(f != k) { dd <- dd + ct[f, g] } }
		}
		bb <- bb + (cc * dd)
	}
	
	myperf <- myperf / (sqrt(aa) * sqrt(bb))
	return(myperf)
}
