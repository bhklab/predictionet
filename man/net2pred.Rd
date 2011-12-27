\name{net2pred}
\alias{net2pred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function fitting a regression model for each gene in the data}
\description{
Function to fit a regression model for each variable in the dataset or alternatively each variable of interest.
}
\usage{
net2pred(net, data, categories, predn, perturbations, method = c("linear", "linear.penalized", "cpt"), seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{net}{ network object. }
	\item{data}{ matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns. }
	\item{categories}{ if this parameter missing, 'data' should be already discretized; otherwise either a single integer or a vector of integers specifying the number of categories used to discretize each variable (data are then discretized using equal-frequency bins) or a list of cutoffs to use to discretize each of the variables in 'data' matrix. If method='bayesnet' and \code{categories} is missing, \code{data} should contain categorical values and the number of categories will determine from the data. }
	\item{predn}{indices or names of variables to fit during network inference. If missing, all the variables will be used for network inference. }
	\item{perturbations}{ matrix of {0,1} specifying whether a gene has been perturbed (e.g., knockdown, overexpression) in some experiments. Dimensions should be the same than \code{data}. }
	\item{method}{ type of predictive model to fit; \code{linear} for linear regression model, \code{linear.penalized} for regularized linear regression model, \code{cpt} for conditional probability tables estimated after discretization of the data. }
	\item{seed}{ set the seed to make the cross-validation and network inference deterministic. }
}
%%\details{
%%
%%}
\value{
	a new network object with the predictive models
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
mynet <- netinf(data=mydata, perturbations=mypert, priors=mypriors, priors.count=TRUE, priors.weight=0.5, maxparents=3, method="regrnet", seed)

net2pred(net=mynet, data=mydata, method="linear")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
