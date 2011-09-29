\name{netinf.cv}
\alias{netinf.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function performing network inference by combining priors and genomic data }
\description{
The function \code{netinf.cv} perform a cross-validation loop and infers a gene network by combining priors and genomic data in each fold. This allows to estimate the predictive ability of the inferred network as well as edge stability.
}
\usage{
netinf.cv(data, categories, perturbations, priors, predn, priors.count = TRUE, priors.weight = 0.5, maxparents = 3, maxparents.push = FALSE, subset, method = c("regrnet", "bayesnet"), ensemble = FALSE, ensemble.maxnsol=3, regrmodel = c("linear", "linear.penalized"), nfold = 10, causal = TRUE, seed, bayesnet.maxcomplexity = 0, bayesnet.maxiter = 100) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns. }
  \item{categories}{ if this parameter missing, 'data' should be already discretize; otherwise either a single integer or a vector of integers specifying the number of categories used to discretize each variable (data are then discretized using equal-frequency bins) or a list of cutoffs to use to discretize each of the variables in 'data' matrix. If method='bayesnet', this parameter should be specified by the user. }
  \item{perturbations}{ matrix of {0, 1} specifying whether a gene has been perturbed (e.g., knockdown, overexpression) in some experiments. Dimensions should be the same than \code{data}. }
  \item{priors}{ matrix of prior information available for gene-gene interaction (parents in rows, children in columns). Values may be probabilities or any other weights (citations count for instance). if priors counts are used the parameter \code{priors.count} should be TRUE so the priors are scaled accordingly. }
  \item{predn}{indices or names of variables to fit during network inference. If missing, all the variables will be used for network inference. }
	\item{priors.count}{ \code{TRUE} if priors specified by the user are number of citations (count) for each interaction, \code{FALSE} if probabilities or any other weight in [0,1] are reported instead. }
  \item{priors.weight}{ real value in [0,1] specifying the weight to put on the priors (0=only the data are used, 1=only the priors are used to infer the topology of the network).}
  \item{maxparents}{ maximum number of parents allowed for each gene. }
 	\item{maxparents.push}{ if TRUE it forces the method to select the maximum number of parents for each variable in the network, FALSE otherwise. For now, this parameters has little effect when method='regrnet' and causal=TRUE. For method='bayesnet', this parameter has an effect the the maximum number of parents is not always chosen, we will fix that in a future release.}
  \item{subset}{ vector of indices to select only subset of the observations. }
  \item{method}{ \code{regrnet} for regression-based network inference, \code{bayesnet} for bayesian network inference with the \code{catnet} package.}
  \item{ensemble}{\code{TRUE} if the ensemble approach should be used, \code{FALSE} otherwise. }
   \item{ensemble.maxnsol}{Number of equivalent solutions chosen at each step.}
	\item{regrmodel }{ type of regression model to fit when \code{regrnet} method is selected. }
	\item{nfold}{ number of folds for the cross-validation. }
	\item{causal}{ 'TRUE' if the causality should be inferred from the data, 'FALSE' otherwise }
	\item{seed}{ set the seed to make the cross-validation and network inference deterministic. }
	\item{bayesnet.maxcomplexity}{ maximum complexity for bayesian network inference, see Details. }
	\item{bayesnet.maxiter}{ maximum number of iterations for bayesian network inference, see Details. }
}
\details{
	\code{bayesnet.maxcomplexity} and \code{bayesnet.maxiter} are parameters to be passed to the network inference method (see \link[catnet]{cnSearchOrder} and \link[catnet]{cnSearchSA} from the \code{catnet} package for more details).
}
\value{
	\item{method }{name of the method used for network inference. }
	\item{topology }{ topology of the model inferred using the entire dataset. }
	\item{topology.coeff }{ if method='regrnet' \code{topology.coeff} contains an adjacency matrix with the coefficients used in the local regression model; parents in rows, children in columns. Additionally the beta_0 values for each model in the first row of the matrix}
	\item{topology.cv }{ topology of the networks inferred at each fold of the cross-validation. }
	\item{topology.coeff.cv }{ if method='regrnet' \code{topology.coeff} contains an adjacency matrix with the coefficients used in the local regression model; parents in rows, children in columns. Additionally the beta_0 values for each model in the first row of the matrix. Inferred at each fold of the cross-validation}
	\item{prediction.score.cv }{ list of prediction scores (R2, NRMSE, MCC) computed at each fold of the cross-validation. }
	\item{edge.stability }{ stability of the edges inferred during cross-validation; only the stability of the edges present in the network inferred using the entire dataset is reported. }
	\item{edge.stability.cv }{stability of the edges inferred during cross-validation. }
	\item{edge.relevance }{ mean relevance score for each across folds in cross-validation. }
	\item{edge.relevance.cv }{ relevance score for each across computed during cross-validation. }
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
## number of fold for cross-validation
res <- netinf.cv(data=mydata, categories=3, perturbations=mypert, priors=mypriors, priors.weight=0.5, method="regrnet", nfold=3, seed=54321)

## MCC for predictions in cross-validation
print(res$prediction.score.cv)

## export network as a 'gml' file that you can import into Cytoscape
rr <- netinf2gml(object=res, file=paste(tempdir(),"predictionet_regrnet",sep=""))

########################
## bayesian network inference
########################
## infer a bayesian network network from data and priors
## number of fold for cross-validation
res <- netinf.cv(data=mydata, categories=3, perturbations=mypert, priors=mypriors, priors.count=TRUE, priors.weight=0.5, method="bayesnet", nfold=3, seed=54321)

## MCC for predictions in cross-validation
print(res$prediction.score.cv)

## export network as a 'gml' file that you can import into Cytoscape
rr <- netinf2gml(object=res, file=paste(tempdir(),"predictionet_bayesnet",sep=""))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
