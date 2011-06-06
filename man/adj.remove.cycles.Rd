\name{adj.remove.cycles}
\alias{adj.remove.cycles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to remove cycles that may be present in a directed graph represented by an adjacency matrix
}
\description{
This function removes cycles that may be present in a directed graph represented by an adjacency matrix, 
}
\usage{
adj.remove.cycles(adjmat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjmat}{
adjacency matrix with positive entries represent evidence for the presence of an edge and entries less or equal than zero represent absence of an edge; parents in row, children in columns.
}
}
\details{
This function may be useful when it comes to generate a bayesian network using a topology identified from an  source of information where cycles are allowed. When cycles are removed, the function tries to keep the most positive entries.
}
\value{
A list of two items
\item{adjmat.acyclic}{ an adjacency matrix without cycles }
\item{adjmat.removed}{ a matrix of booleans representing the edges that have been removed from the original adjacency matrix to make it acyclic }
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%%}
\examples{
set.seed(54321)
xx <- matrix(sample(c(0,1), 100, replace=TRUE), nrow=10, ncol=10)
adj.remove.cycles(adjmat=xx)
}
\keyword{ graphs }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
