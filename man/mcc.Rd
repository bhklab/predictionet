\name{mcc}
\alias{mcc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the Matthews Correlation Coefficient (MCC) in a classification framework
}
\description{
This function computes the Matthews Correlation Coefficient (MCC) in a classification framework.
}
\usage{
mcc(ct, nbcat = nrow(ct))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ct}{
contingency table
}
  \item{nbcat}{
number of categories
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
MCC estimate
}
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
\keyword{ classif }% __ONLY ONE__ keyword per line
