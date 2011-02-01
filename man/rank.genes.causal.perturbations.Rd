\name{rank.genes.causal.perturbations}
\alias{rank.genes.causal.perturbations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function ranking the genes based on mrmr+prior+causality for each gene as a target
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rank.genes.causal.perturbations(priors, data, perturbations, predn, priors.weight, maxparents, causal = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{priors}{
%%     ~~Describe \code{priors} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{perturbations}{
%%     ~~Describe \code{perturbations} here~~
}
	\item{predn}{indices variables to fit during network inference. If missing, all the variables will be used for network inference. }
  \item{priors.weight}{
%%     ~~Describe \code{priors.weight} here~~
}
  \item{maxparents}{
%%     ~~Describe \code{maxparents} here~~
}
  \item{causal}{
%%     ~~Describe \code{maxparents} here~~
}
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
