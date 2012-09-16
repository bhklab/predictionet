\name{predictionet.press.statistic}
\alias{predictionet.press.statistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function computing the press statistic for all target variables in topology}
\description{
The function \code{predictionet.press.statistic} computes the press statistic for all target variables in the provided topology.
}
\usage{
predictionet.press.statistic(topo,data,ensemble=FALSE,perturbations=NULL) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{topo}{ adjacency matrix of {0,1} indicating whether two variables are connected}
  \item{data}{ matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns. }
  \item{perturbations}{ matrix of {0, 1} specifying whether a gene has been perturbed (e.g., knockdown, overexpression) in some experiments. Dimensions should be the same than \code{data}. }
  \item{ensemble}{ \code{TRUE} if the ensemble approach should be used, \code{FALSE} otherwise. }
}
%%\details{
%%
%%}
\value{
	A vector of press statistics, one for every target variable.
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
