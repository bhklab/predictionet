\name{predictionet.stability.cv}
\alias{predictionet.stability.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function inferring networks in cross-validation }
\description{
The function \code{predictionet.stability.cv} infers networks in cross-validation (compared to netinf.cv no regression is carried out, thus less computational cost but no prediction scores)

}
\usage{
predictionet.stability.cv(data, categories, perturbations, priors, predn, priors.count = TRUE, priors.weight = 0.5, maxparents = 3, subset, method = c("regrnet", "bayesnet"), ensemble = FALSE, ensemble.maxnsol = 3, nfold = 10, causal = TRUE, seed, bayesnet.maxcomplexity = 0, bayesnet.maxiter = 100) 
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
  \item{subset}{ vector of indices to select only subset of the observations. }
  \item{method}{ \code{regrnet} for regression-based network inference, \code{bayesnet} for bayesian network inference with the \code{catnet} package.}
  \item{ensemble}{\code{TRUE} if the ensemble approach should be used, \code{FALSE} otherwise. }
   \item{ensemble.maxnsol}{Number of equivalent solutions chosen at each step.}
	\item{nfold}{ number of folds for the cross-validation. }
	\item{causal}{ 'TRUE' if the causality should be inferred from the data, 'FALSE' otherwise }
	\item{seed}{ set the seed to make the cross-validation and network inference deterministic. }
	\item{bayesnet.maxcomplexity}{ maximum complexity for bayesian network inference, see Details. }
	\item{bayesnet.maxiter}{ maximum number of iterations for bayesian network inference, see Details. }}
%%\details{
%%
%%}
\value{
	\item{method }{name of the method used for network inference. }
	\item{topology }{ topology of the model inferred using the entire dataset. }
	\item{topology.cv }{ topology of the networks inferred at each fold of the cross-validation. }
	\item{edge.stability }{ stability of the edges inferred during cross-validation; only the stability of the edges present in the network inferred using the entire dataset is reported. }
	\item{edge.stability.cv }{stability of the edges inferred during cross-validation. }
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
res <-netinf.cv(data=mydata, categories=3, perturbations=mypert, priors=mypriors, priors.weight=0.5, method="regrnet", nfold=3, seed=54321)

## MCC for predictions in cross-validation
print(res$prediction.score.cv)

## export network as a 'gml' file that you can import into Cytoscape
\dontrun{rr <- netinf2gml(object=res, file="predictionet_regrnet")}

########################
## bayesian network inference
########################
## infer a bayesian network network from data and priors
## number of fold for cross-validation
\dontrun{res <- netinf.cv(data=mydata, categories=3, perturbations=mypert, priors=mypriors, priors.count=TRUE, priors.weight=0.5, method="bayesnet", nfold=3, seed=54321)}

## MCC for predictions in cross-validation
\dontrun{print(res$prediction.score.cv)}

## export network as a 'gml' file that you can import into Cytoscape
\dontrun{rr <- netinf2gml(object=res, file="predictionet_bayesnet")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
