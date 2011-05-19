\name{fit.regrnet.causal}
\alias{fit.regrnet.causal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function fitting a regression model for each gene in the data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fit.regrnet.causal(data, perturbations, priors, predn, maxparents = 3, priors.weight = 0.5, regrmodel = c("linear", "linear.penalized"), causal = TRUE, seed = 54321)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{perturbations}{
%%     ~~Describe \code{perturbations} here~~
}
  \item{priors}{
%%     ~~Describe \code{priors} here~~
}
	\item{predn}{indices or names of variables to fit during network inference. If missing, all the variables will be used for network inference. }
  \item{maxparents}{
%%     ~~Describe \code{maxparents} here~~
}
  \item{priors.weight}{
%%     ~~Describe \code{priors.weight} here~~
}
	\item{regrmodel }{ type of regression model to fit when \code{regrnet} method is selected. }
	\item{causal}{ 'TRUE' if the causality should be inferred from the data, 'FALSE' otherwise }
	\item{seed}{ set the seed to make the network inference deterministic. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
