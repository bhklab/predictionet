\name{build.regression.regrnet}
\alias{build.regression.regrnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function fitting a regression model for each gene in the data}
\description{
The function \code{build.regression.regrnet}  fits a regression model for each variable in the dataset or alternatively each variable specified in predn.
}
\usage{
build.regression.regrnet(data, predn, topo , regrmodel=c("linear", "linear.penalized"),ensemble=FALSE, seed=54321)  
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns. }
    \item{predn}{indices or names of variables to fit during network inference. If missing, all the variables will be used for network inference. }
	  \item{topo}{ adjacency matrix of {0,1} indicating whether two variables are connected}
  	\item{regrmodel }{ type of regression model to fit when \code{regrnet} method is selected. }
  \item{ensemble}{ \code{TRUE} if the ensemble approach should be used, \code{FALSE} otherwise. }
	\item{seed}{ set the seed to make the cross-validation and network inference deterministic. }
}
%%\details{
%%
%%}
\value{
	A list of regression models.
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
