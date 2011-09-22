\name{eval.network}
\alias{eval.network}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function computing the f1-score, comparing an inferred topology with a given topology
}
\description{
This function computes the f1-score for an inferred topology using a topology provided by the user.
}
\usage{
eval.network(net,true.net)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
Inferred topology, an edge between to variables X and Y corresponds to net[X,Y]=1.}
  \item{true.net}{
topology the user wants to compare the inferred topology with, e.g. the true network using generated datasets. An edge between to variables X and Y corresponds to net[X,Y]=1.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
The computed f1-score, defined as 2*TP/(2*TP+FN+FP)
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
## regression-based network inference
res <- netinf(data=mydata, categories=3, priors=mypriors, priors.weight=0.5, method="regrnet", seed=54321)

## extract adjacency matrix from inferred network
net2topo(net=res)
## with coefficients
net2topo(net=res, coefficients=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
