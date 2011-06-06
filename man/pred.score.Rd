\name{pred.score}
\alias{pred.score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function computing performance of prediction; methods include r2, nrmse and mcc
}
\description{
This function computes prediction performance; methods include r2, nrmse and mcc.
}
\usage{
pred.score(data, pred, categories, method = c("r2", "nrmse", "mcc"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{pred}{
%%     ~~Describe \code{pred} here~~
}
  \item{categories}{ if this parameter missing, 'data' should be already discretize; otherwise either a single integer or a vector of integers specifying the number of categories used to discretize each variable (data are then discretized using equal-frequency bins) or a list of cutoffs to use to discretize each of the variables in 'data' matrix. If method='bayesnet', this parameter should be specified by the user. }
  \item{method}{
%%     ~~Describe \code{method} here~~
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
A vector of performance scores, one for each node
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

\seealso{
\link[predictionet]{netinf.predict}
}
\examples{
set.seed(54321)
xx <- runif(100)
## R2
pred.score(data=xx, pred=xx+rnorm(100)/10, method="r2")
## NRMSE
pred.score(data=xx, pred=xx+rnorm(100)/10, method="nrmse")
## MCC
pred.score(data=xx, pred=xx+rnorm(100)/10, categories=3, method="mcc")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ classif }
\keyword{ regression }