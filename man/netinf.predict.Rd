\name{netinf.predict}
\alias{netinf.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to make prediction of a node values given its parents using an inferred network
}
\description{
This function predict the value of a node given its parents using an inferred network
}
\usage{
netinf.predict(net, data, categories, perturbations, subset, predn, ensemble=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
%%     ~~Describe \code{net} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{categories}{
%%     ~~Describe \code{categories} here~~
}
  \item{perturbations}{
%%     ~~Describe \code{perturbations} here~~
}
  \item{subset}{
%%     ~~Describe \code{subset} here~~
}
	\item{predn}{indices or names of variables to fit during network inference. If missing, all the variables will be used for network inference. }
		  \item{ensemble}{\code{TRUE} if the ensemble approach should be used, \code{FALSE} otherwise. }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
matrix of predicted values
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Benjamin Haibe-Kains, Catharina Olsen
}
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
## number of genes to select for the analysis
genen <- 10
## select only the top genes
goi <- dimnames(annot.ras)[[1]][order(abs(log2(annot.ras[ ,"fold.change"])), decreasing=TRUE)[1:genen]]
mydata <- data.ras[ , goi, drop=FALSE]
myannot <- annot.ras[goi, , drop=FALSE]
mypriors <- priors.ras[goi, goi, drop=FALSE]
mydemo <- demo.ras
## infer global network from data and priors
mynet <- netinf(data=mydata, priors=mypriors, priors.count=TRUE, priors.weight=0.5, maxparents=3, method="regrnet", seed=54321)

## predict gene expression of the first gene
mypreds <- netinf.predict(net=mynet, data=mydata, predn=goi[1])[ ,goi[1]]
## root mean squared error (RMSE)
sqrt(mean((mydata[ ,goi[1]] - mypreds)^2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classif }
\keyword{ regression }% __ONLY ONE__ keyword per line
