## function to make the network model predictive either by fitting local regression models or by estimating CPTs
`net2pred` <- 
function(net, data, categories, predn, perturbations, method=c("linear", "linear.penalized", "cpt"), seed) {
	method <- match.arg(method)
	if(missing(perturbations) || is.null(perturbations)) {
		perturbations <- matrix(FALSE, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	} else {
		if(nrow(perturbations) == 1) {
			perturbations[1, ] <- as.logical(perturbations[1, ])
		} else { perturbations <- apply(perturbations, 2, as.logical) }
		dimnames(perturbations) <- dimnames(data)
	}
	if(length(method) > 1) { stop("only one prediction can be specified!") }
	switch(method, 
	"linear"=, 
	"linear.penalized"={
		return(.build.regression.regrnet(net=net, data=data, predn=predn, perturbations=perturbations, regrmodel=method, seed=seed))
	}, 
	"cpt"={
		stop("cpt method is not implemented yet!")
	}, 
	stop("no default parameter for method!")
	)
}
