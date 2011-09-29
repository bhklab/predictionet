\name{netinf}
\alias{netinf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function performing network inference by combining priors and genomic data }
\description{
Main function of the \code{predictionet} package, \code{netinf} infers a gene network by combining priors and genomic data. The two main network inference methodologies implemented so far are the bayesian and regression-based inferences.
}
\usage{
netinf(data, categories, perturbations, priors, predn, priors.count = TRUE, priors.weight = 0.5, maxparents = 3, maxparents.push = FALSE, subset, method = c("regrnet", "bayesnet"), ensemble = FALSE, regrmodel = c("linear", "linear.penalized"), ensemble.model = c("full","best"), ensemble.maxnsol = 3, causal=TRUE, seed=54321, bayesnet.maxcomplexity=0, bayesnet.maxiter=100) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns. }
  \item{categories}{ if this parameter missing, 'data' should be already discretized; otherwise either a single integer or a vector of integers specifying the number of categories used to discretize each variable (data are then discretized using equal-frequency bins) or a list of cutoffs to use to discretize each of the variables in 'data' matrix. If method='bayesnet' and \code{categories} is missing, \code{data} should contain categorical values and the number of categories will determine from the data. }
  \item{perturbations}{ matrix of {0,1} specifying whether a gene has been perturbed (e.g., knockdown, overexpression) in some experiments. Dimensions should be the same than \code{data}. }
  \item{priors}{ matrix of prior information available for gene-gene interaction (parents in rows, children in columns). Values may be probabilities or any other weights (citations count for instance). if priors counts are used the parameter \code{priors.count} should be TRUE so the priors are scaled accordingly. }
	\item{predn}{indices or names of variables to fit during network inference. If missing, all the variables will be used for network inference. Note that for bayesian network inference (method='bayesnet') this parameter is ignored and a network will be generated using all the variables.}
	\item{priors.count}{ \code{TRUE} if priors specified by the user are number of citations (count) for each interaction, \code{FALSE} if probabilities or any other weight in [0,1] are reported instead. }
  \item{priors.weight}{ real value in [0,1] specifying the weight to put on the priors (0=only the data are used, 1=only the priors are used to infer the topology of the network).}
  \item{maxparents}{ maximum number of parents allowed for each gene. }
 	\item{maxparents.push}{ if TRUE it forces the method to select the maximum number of parents for each variable in the network, FALSE otherwise. For now, this parameters has little effect when method='regrnet' and causal=TRUE. For method='bayesnet', this parameter has an effect the the maximum number of parents is not always chosen, we will fix that in a future release.}
 \item{subset}{ vector of indices to select only subset of the observations. }
  \item{method}{ \code{regrnet} for regression-based network inference, \code{bayesnet} for bayesian network inference with the \code{catnet} package.}
  \item{ensemble}{ \code{TRUE} if the ensemble approach should be used, \code{FALSE} otherwise. }
	\item{regrmodel }{ type of regression model to fit when \code{regrnet} method is selected. }
	\item{ensemble.model}{ Could be either \code{full} or \code{best} depending how the equivalent networks are selected to be included in the ensemble network: for \code{full} bootstrapping is used to identify all the statistically equivalent networks, it \code{best} only the top \code{ensemble.maxnsol} are considered at each step of the feature selection.}
	\item{ensemble.maxnsol}{ maximum number of solutions to consider at each step of the feature selection for the method=\code{ensemble.regrnet}, default is 3.}
	\item{causal}{ 'TRUE' if the causality should be inferred from the data, 'FALSE' otherwise }
	\item{seed}{ set the seed to make the network inference deterministic. }
	\item{bayesnet.maxcomplexity}{ maximum complexity for bayesian network inference, see Details. }
	\item{bayesnet.maxiter}{ maximum number of iterations for bayesian network inference, see Details. }
}
\details{
\code{bayesnet.maxcomplexity} and \code{bayesnet.maxiter} are parameters to be passed to the network inference method (see \link[catnet]{cnSearchOrder} and \link[catnet]{cnSearchSA} from the \code{catnet} package for more details).

Relevance score is either MRMR scores if causal=FALSE or causality score if causal=FALSE.
}
\value{
	\item{method }{name of the method used for network inference. }
	\item{topology }{ adjacency matrix representing the topology of the inferred network; parents in rows, children in columns. }
	\item{topology.coeff }{ if method='regrnet' \code{topology.coeff} contains an adjacency matrix with the coefficients used in the local regression model; parents in rows, children in columns. Additionally the beta_0 values for each model in the first row of the matrix}
	\item{edge.relevance }{ relevance score for each edge (see Details). }
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{ Benjamin Haibe-Kains, Catharina Olsen }
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
## load gene expression data for colon cancer data, list of genes related to RAS signaling pathway and the corresponding priors
data(expO.colon.ras)
## create matrix of perturbations (no perturbations in this dataset)
pert <- matrix(0, nrow=nrow(data.ras), ncol=ncol(data.ras), dimnames=dimnames(data.ras))

## number of genes to select for the analysis
genen <- 10
## select only the top genes
goi <- dimnames(annot.ras)[[1]][order(abs(log2(annot.ras[ ,"fold.change"])), decreasing=TRUE)[1:genen]]
mydata <- data.ras[ , goi, drop=FALSE]
myannot <- annot.ras[goi, , drop=FALSE]
mypriors <- priors.ras[goi, goi, drop=FALSE]
mydemo <- demo.ras
mypert <- pert[ , goi, drop=FALSE]

########################
## regression-based network inference
########################
## infer global network from data and priors
mynet <- netinf(data=mydata, perturbations=mypert, priors=mypriors, priors.count=TRUE, priors.weight=0.5, maxparents=3, method="regrnet", seed=54321)

## plot network topology
mytopo <- mynet$topology
library(network)
xnet <- network(x=mytopo, matrix.type="adjacency", directed=TRUE, loops=FALSE, vertex.attrnames=dimnames(mytopo)[[1]])
plot.network(x=xnet, displayisolates=TRUE, displaylabels=TRUE, boxed.labels=FALSE, label.pos=0, arrowhead.cex=2, vertex.cex=4, vertex.col="royalblue", jitter=FALSE, pad=0.5)

## export network as a 'gml' file that you can import into Cytoscape
rr <- netinf2gml(object=mynet, file=paste(tempdir(),"/predictionet_regrnet",sep=""))

## predictions
mypreds <- netinf.predict(net=mynet, data=mydata, perturbations=mypert)
print(pred.score(data=mydata, pred=mypreds, categories=3, method="mcc"))

########################
## bayesian network inference
########################
## discretize gene expression values in three categories
categories <- rep(3, ncol(mydata))
## estimate the cutoffs (tertiles) for each gene
cuts.discr <- lapply(apply(rbind("nbcat"=categories, mydata), 2, function(x) { y <- x[1]; x <- x[-1]; return(list(quantile(x=x, probs=seq(0, 1, length.out=y+1), na.rm=TRUE)[-c(1, y+1)])) }), function(x) { return(x[[1]]) })
mydata <- data.discretize(data=mydata, cuts=cuts.discr)

## infer a bayesian network network from data and priors
mynet <- netinf(data=mydata, perturbations=mypert, priors=mypriors, priors.count=TRUE, priors.weight=0.5, maxparents=3, method="bayesnet", seed=54321)

## plot network topology
mytopo <- net2topo(net=mynet)
library(network)
xnet <- network(x=mytopo, matrix.type="adjacency", directed=TRUE, loops=FALSE, vertex.attrnames=dimnames(mytopo)[[1]])
plot.network(x=xnet, displayisolates=TRUE, displaylabels=TRUE, boxed.labels=FALSE, label.pos=0, arrowhead.cex=2, vertex.cex=4, vertex.col="royalblue", jitter=FALSE, pad=0.5)

## export network as a 'gml' file that you can import into Cytoscape
rr <- netinf2gml(object=mynet, file=paste(tempdir(),"/predictionet_bayesnet",sep=""))

## predictions
mypreds <- netinf.predict(net=mynet, data=mydata, perturbations=mypert)
print(pred.score(data=mydata, pred=mypreds, method="mcc"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
