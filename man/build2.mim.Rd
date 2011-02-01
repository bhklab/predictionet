\name{build2.mim}
\alias{build2.mim}

\title{Build Mutual Information Matrix (R implementation)}
\description{
	  \code{build2.mim} takes the dataset as input and computes the 
	  mutual information between all pair of variables according 
	  to the mutual inforamtion estimator \code{estimator}. 
	  The results are saved in the mutual information matrix (MIM), a square 
      matrix whose (i,j) element is the mutual information between variables
      \eqn{X_i}{Xi} and \eqn{X_j}{Xj}.
}
\usage{
build2.mim(dataset, estimator = c("pearson", "spearman", "kendall"))
}
\arguments{
  \item{dataset}{data.frame containing gene expression data or any dataset where columns contain variables/features and rows contain outcomes/samples.}
  \item{estimator}{The name of the entropy estimator to be used. The package implements four estimators for discrete data: three estimators based on correlation (can be used with continuous data): "pearson","spearman","kendall"(default:"spearman") - see details.
	}

}
\value{ \code{build2.mim} returns the mutual information matrix.
}

\details{
This R code was extracted from the \code{minet} package developed and maintaing by Patrick Meyer (see reference).

Estimators:
\itemize{
	  \item "pearson" :  This computes mutual information for normally distributed variable.
	  \item "spearman" : This computes mutual information for normally distributed variable using Spearman's correlation instead of Pearson's correlation.
	  \item "kendall" : This computes mutual information for normally distributed variable using Kendall's correlation instead of Pearson's correlation.
}
}
\author{
  Patrick E. Meyer, Frederic Lafitte, Gianluca Bontempi
}
\references{
  Patrick E. Meyer, Frederic Lafitte, and Gianluca Bontempi. 
 minet: A R/Bioconductor Package for Inferring Large Transcriptional Networks Using Mutual Information. BMC Bioinformatics, Vol 9,
 2008.
 
 J. Beirlant, E. J. Dudewica, L. Gyofi, and E. van der Meulen. Nonparametric 
 entropy estimation : An overview. Journal of Statistics, 1997.

 Jean Hausser. Improving entropy estimation and the inference of genetic regulatory networks.
 Master thesis of the National Institute of Applied Sciences of Lyon, 2006.
}
\seealso{\code{\link[minet]{clr}}, \code{\link[minet]{aracne}}, \code{\link[minet]{mrnet}}, \code{\link[minet]{mrnetb}}}
\examples{
## load gene expression data for colon cancer data, list of genes related to RAS signaling pathway and the corresponding priors
data(expO.colon.ras)
## let's compute the MIM
mim <- build2.mim(data.ras, estimator="pearson")
}
\keyword{misc}
