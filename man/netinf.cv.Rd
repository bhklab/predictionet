\name{netinf.cv}
\alias{netinf.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function performing network inference by combining priors and genomic data }
\description{
The function \code{netinf.cv} perform a cross-validation loop and infers a gene network by combining priors and genomic data in each fold. This allows to estimate the predictive ability of the inferred network as well as edge stability.
}
\usage{
netinf.cv(data, categories, perturbations, priors, predn, priors.count = TRUE, priors.weight = 0.5, maxparents = 3, subset, method = c("regrnet", "regrnet.ensemble", "bayesnet"), nfold=10, causal = TRUE, seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns. }
  \item{categories}{ if this parameter missing, 'data' should be already discretize; otherwise either a single integer or a vector of integers specifying the number of categories used to discretize each variable (data are then discretized using equal-frequency bins) or a list of cutoffs to use to discretize each of the variables in 'data' matrix. If method='bayesnet', this parameter should be specified by the user. }
  \item{perturbations}{ matrix of {0, 1} specifying whether a gene has been perturbed (e.g., knockdown, overexpression) in some experiments. Dimensions should be the same than \code{data}. }
  \item{priors}{ matrix of prior information available for gene-gene interaction (parents in rows, children in columns). Values may be probabilities or any other weights (citations count for instance). }
  \item{predn}{indices or names of variables to fit during network inference. If missing, all the variables will be used for network inference. }
	\item{priors.count}{ \code{TRUE} if priors specified by the user are number of citations (count) for each interaction, \code{FALSE} if probabilities are reported instead. }
  \item{priors.weight}{ real value in [0, 1] specifying the weight to put on the priors (0=only the data are used, 1=only the priors are used to infer the topology of the network). If \code{priors.weight} is missing it will be optimized gene by gene in an automatic way. }
  \item{maxparents}{ maximum number of parents allowed for each gene. }
  \item{subset}{ vector of indices to select only subset of the observations. }
  \item{method}{ \code{bayesnet} for bayesian network inference with the \code{catnet} package (not implemented yet), \code{regrnet} for regression-based network inference, \code{regrnet.ensemble} for ensemble regression-based network inference. }
	\item{nfold}{ number of folds for the cross-validation. }
	\item{causal}{ 'TRUE' if the causality should be inferred from the data, 'FALSE' otherwise }
	\item{seed}{ set the seed to make the cross-validation and network inference deterministic. }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
	\item{net }{ model object of the network inferred using the entire dataset. }
	\item{net.cv }{ list of network models fitted at each fold of the cross-validation. }
	\item{topology }{ topology of the model inferred using the entire dataset. }
	\item{topology.cv }{ topology of the networks inferred at each fold of the cross-validation. }
	\item{prediction.score.cv }{ list of prediction scores (R2, NRMSE, MCC) computed at each fold of the cross-validation. }
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

## number of fold for cross-validation
res <- netinf.cv(data=mydata, categories=3, perturbations=mypert, priors=mypriors, method="regrnet", nfold=3, seed=54321)

## MCC for predictions in cross-validation
print(res$prediction.score.cv)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
