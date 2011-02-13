\name{regrnet2topo}
\alias{regrnet2topo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function transforming a regression model into an adjacency matrix: parents in rows, children in columns.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
regrnet2topo(net, coefficients = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
	network model extracted from the 'net' item of a network model object.
}
  \item{coefficients}{
	logical specifying whether the adjacency matrix contains the coefficients fitted in each local regression model or 0/1. 'TRUE' for the coefficients, 'FALSE' for generating a matrix of binary indicators (0 = no edge, 1 = edge linking gene i to gene j).
}
}
\details{
Note that a '0' in the binary adjacency matrix represents the absence of an edge BEFORE the fitting of the local (penalized) regression model while in the matrix containing the coefficients '0' means that this gene expression is not used in the local regression model, either because it was not selected for regression or because the coefficient is 0 in the fitted regression model.
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
