\name{predictionet-package}
\alias{predictionet-package}
\alias{predictionet}
\docType{package}
\title{
Inference for predictive networks designed for (but not limited to) genomic data
}
\description{
This package contains a set of functions related to network inference combining genomic data and prior information extracted from biomedical literature and structured biological databases.The main function is able to generate networks using bayesian or regression-based inference methods; while the former is limited to < 100 of variables, the latter may infer network with hundreds of variables. Several statistics at the edge and node levels have been implemented (edge stability, predictive ability of each node, ...) in order to help the user to focus on high quality subnetworks. Ultimately, this package is used in the 'Predictive Networks' web application developed by the Dana-Farber Cancer Institute in collaboration with
}
\details{
\tabular{ll}{
Package: \tab predictionet\cr
Type: \tab Package\cr
Version: \tab \Sexpr{packageDescription("predictionet")$Version}\cr
Date: \tab \Sexpr{packageDescription("predictionet")$Date}\cr
License: \tab Artistic 2.0\cr
LazyLoad: \tab yes\cr
}
%%~~ An overview of how to use the package, including the most important ~~
%%~~ functions ~~
}
\author{
Benjamin Haibe-Kains, Catharina Olsen, Gianluca Bontempi, John Quackenbush

- Computational Biology and Functional Genomics, Dana-Farber Cancer Institute, Boston, MA, USA

\url{http://compbio.dfci.harvard.edu/}

- Center for Cancer Computational Biology, Dana-Farber Cancer Institute, Boston, MA, USA

\url{http://cccb.dfci.harvard.edu/index.html}

- Machine Learning Group (MLG), Universite Libre de Bruxelles, Bruxelles, Belgium 

\url{http://www.ulb.ac.be/di/mlg/}


\bold{Maintainer}:
Benjamin Haibe-Kains

\email{bhaibeka@jimmy.harvard.edu}

\email{bhaibeka@ulb.ac.be}

Catharina Olsen

\email{colsen@ulb.ac.be}

}
%%\references{
%%~~ Literature or other references for background information ~~
%%}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%%~~ the R documentation directory ~~
\keyword{ package }
%%\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%%}
%%\examples{
%%~~ simple examples of the most important functions ~~
%%}
