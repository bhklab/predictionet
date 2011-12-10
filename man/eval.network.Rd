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
eval.network(topo, true.topo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{topo}{
Inferred topology, an edge between to variables X and Y corresponds to net[X,Y]=1.}
  \item{true.topo}{
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
%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
