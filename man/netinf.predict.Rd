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
netinf.predict(net, data, categories, perturbations, subset, predn, method=c("linear", "linear.penalized", "cpt"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{ a network object with local regression models. }
  \item{data}{ matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns. }
  \item{categories}{ if this parameter missing, 'data' should be already discretize; otherwise either a single integer or a vector of integers specifying the number of categories used to discretize each variable (data are then discretized using equal-frequency bins) or a list of cutoffs to use to discretize each of the variables in 'data' matrix. If method='bayesnet', this parameter should be specified by the user. }
  \item{perturbations}{ matrix of {0, 1} specifying whether a gene has been perturbed (e.g., knockdown, overexpression) in some experiments. Dimensions should be the same than \code{data}. }
 \item{subset}{ vector of indices to select only subset of the observations. }
	\item{predn}{indices or names of variables to fit during network inference. If missing, all the variables will be used for network inference. }
	\item{method}{ \code{regrnet} for regression-based network inference, \code{bayesnet} for bayesian network inference with the \code{catnet} package.}
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
mynet <- net2pred(net=mynet, data=mydata, method="linear")

## predict gene expression of the first gene
mypreds <- netinf.predict(net=mynet, data=mydata, predn=goi[1])[ ,goi[1]]
## root mean squared error (RMSE)
nrmse <- sqrt(mean((mydata[ ,goi[1]] - mypreds)^2))
## R2
r2 <- cor(mydata[ ,goi[1]], mypreds)^2
plot(mydata[ ,goi[1]], mypreds, xlab="Observed gene expression", ylab="Predicted gene expression")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classif }
\keyword{ regression }% __ONLY ONE__ keyword per line
