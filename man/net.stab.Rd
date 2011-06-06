\name{net.stab}
\alias{net.stab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function evaluating the stability of a network inferred from data with a specific method
}
\description{
This function evaluates the stability of a network inferred from data with a specific method.
}
\usage{
net.stab(data, categories, perturbations, priors, priors.count = TRUE, priors.weight = 0.5, maxparents = 3, subset, method = c("bayesnet", "regrnet"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{categories}{
%%     ~~Describe \code{categories} here~~
}
  \item{perturbations}{
%%     ~~Describe \code{perturbations} here~~
}
  \item{priors}{
%%     ~~Describe \code{priors} here~~
}
  \item{priors.count}{
%%     ~~Describe \code{priors.count} here~~
}
  \item{priors.weight}{
%%     ~~Describe \code{priors.weight} here~~
}
  \item{maxparents}{
%%     ~~Describe \code{maxparents} here~~
}
  \item{subset}{
%%     ~~Describe \code{subset} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
%%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%%}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ graphs }% __ONLY ONE__ keyword per line
