\name{adj.descent}
\alias{adj.descent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to identify all children of a parent
}
\description{
This function uses a depth-first search algorithm to identify all the children (and their corresponding depth) of a node.
}
\usage{
adj.descent(adjmat, from, depth)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjmat}{
adjacency matrix; parents in rows, children in columns
}
  \item{from}{
name or index of the parent node
}
  \item{depth}{
numeric representing the depth of the search for the children; if \code{depth} is missing or negative, all children are returned
}
}
\details{
The algorithm is based on the depth-first search.
}
\value{
two-column matrix containing the names of the children in the first column and their corresponding depth in the descent in the second column
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
\examples{
## check whether a list of two nodes are children of another node
set.seed(54321)
mytopo <- matrix(sample(0:1, 100, replace=TRUE, prob=c(0.7,0.3)), nrow=10, dimnames=list(LETTERS[1:10], LETTERS[1:10]))
adj.descent(adjmat=mytopo, from="A")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
