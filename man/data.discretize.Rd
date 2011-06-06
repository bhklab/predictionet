\name{data.discretize}
\alias{data.discretize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to discretize data based on user specified cutoffs
}
\description{
This function enable  discretization of data based on cutoffs specified by the users
}
\usage{
data.discretize(data, cuts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
matrix of continuous or categorical values (gene expressions for example); observations in rows, features in columns.
}
  \item{cuts}{
list of cutoffs for each variable.
}
}
\details{
This function is discretizing the continuous value in \code{data} using the cutoffs specified in \code{cuts} to create categories represented by increasing integers in {1,2,...n} where n is the maximum number of categories in the dataset.
}
\value{
a matrix of categorical values where categories are \{1,2,..,n\} depending on the list of cutoffs specified in \code{cuts}; observations in rows, features in columns.
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

\seealso{
\link[infotheo]{discretize}
}
\examples{
## load gene expression data for colon cancer data, list of genes related to RAS signaling pathway and the corresponding priors
data(expO.colon.ras)
## discretize the data in 3 categories
categories <- rep(3, ncol(data.ras))
## estimate the cutoffs (tertiles) for each gene
cuts.discr <- lapply(apply(rbind("nbcat"=categories, data.ras), 2, function(x) { y <- x[1]; x <- x[-1]; return(list(quantile(x=x, probs=seq(0, 1, length.out=y+1), na.rm=TRUE)[-c(1, y+1)])) }), function(x) { return(x[[1]]) })
data.ras.bin <- data.discretize(data=data.ras, cuts=cuts.discr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
