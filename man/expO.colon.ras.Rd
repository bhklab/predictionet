\name{expO.colon.ras}
\alias{expO.colon.ras}
\alias{demo.ras}
\alias{data.ras}
\alias{annot.ras}
\alias{priors.ras}
\docType{data}
\title{
Gene expression, annotations, clinical data and priors for the colon cancer tumors collected by the expression project for oncology (expO).
}
\description{
This dataset contains (part of) the gene expression, annotations and clinical data as published by the expO project (\url{http://www.intgen.org/expo/}). Genes related to KRAS mutations were retrieved from Bild et al, Nature, 2006. Only genes with known gene symbols were selected resulting in a dataset of 292 human colon tumors and 259 RAS-related genes.
}
\usage{data(expO.colon.ras)}
\format{
\code{\link[predictionet]{expO.colon.ras}} is a dataset containing four matrices:
 \describe{
	\item{demo.ras}{clinical information of the colon cancer patients whose tumors were hybridized.}
    \item{data.ras}{matrix containing expression of genes related to RAS.}
    \item{annot.ras}{matrix containing annotations of the genes related to RAS.}
	\item{priors.ras}{matrix of priors counts for all the genes related to RAS. Each value represents the number of times an interaction was observed for a specific pair of genes (parents in rows, children in columns).}
}
}
\details{
The microarray platform used in the expO project is the Affymetrix HG-U133PLUS2 GeneChip.
}
\source{
\url{https://expo.intgen.org/geo/}

\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2109}
}
\references{
\url{http://www.intgen.org/expo/}

Bild AH, Yao G, Chang JT, Wang Q, Potti A, Chasse D, Joshi MB, Harpole D, Lancaster JM, Berchuck A, Olson JA Jr, Marks JR, Dressman HK, West M, Nevins JR. (2006) "Oncogenic pathway signatures in human cancers as a guide to targeted therapies", \emph{Nature}, \bold{439}(7074):274-275.
}
\examples{
data(expO.colon.ras)
}
\keyword{data}
