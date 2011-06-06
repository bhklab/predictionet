\name{net2topo}
\alias{net2topo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function transforming a network model into an adjacency matrix: parents in rows, children in columns
}
\description{
This function transforms a network model into an adjacency matrix: parents in rows, children in columns.
}
\usage{
net2topo(net, coefficients=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
%%     ~~Describe \code{net} here~~
}
  \item{coefficients}{
%%     ~~Describe \code{coefficients} here~~
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Matrix of binary values (0=absence of an edge and 1=presence of an edge) or regression coefficients as estimated in the local regression models (input=parents nodes, output=target/child node)
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
