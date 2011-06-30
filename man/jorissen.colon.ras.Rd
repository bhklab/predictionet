\name{jorissen.colon.ras}
\alias{jorissen.colon.ras}
\alias{demo2.ras}
\alias{data2.ras}
\alias{annot2.ras}
\alias{priors2.ras}
\docType{data}
\title{
Gene expression, annotations, clinical data and priors for the colon cancer tumors collected by Jorissen and colleagues in 2009.
}
\description{
This dataset contains (part of) the gene expression, annotations and clinical data as published by Jorissen and colleagues in 2009. Genes related to KRAS mutations were retrieved from Bild et al, Nature, 2006. Only genes with known gene symbols were selected resulting in a dataset of 290 human colon tumors and 259 RAS-related genes.
}
\usage{data(jorissen.colon.ras)}
\format{
\code{\link[predictionet]{jorissen.colon.ras}} is a dataset containing four matrices:
\describe{
\item{demo2.ras}{clinical information of the colon cancer patients whose tumors were hybridized.}
\item{data2.ras}{matrix containing expression of genes related to RAS.}
\item{annot2.ras}{matrix containing annotations of the genes related to RAS.}
\item{priors2.ras}{matrix of priors counts for all the genes related to RAS. Each value represents the number of times an interaction was observed for a specific pair of genes (parents in rows, children in columns).}
}
}
\details{
The microarray platform used in Jorissen's dataset is the Affymetrix HG-U133PLUS2 GeneChip.
}
\source{
\url{http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE14333}
}
\references{
Jorissen RN, Gibbs P, Christie M, Prakash S, Lipton L, Desai J, Kerr D, Aaltonen LA, Arango D, Kruhoffer M, Orntoft TF, Andersen CL, Gruidl M, Kamath VP, Eschrich S, Yeatman TJ, Sieber OM. (2009) "Metastasis-Associated Gene Expression Changes Predict Poor Outcomes in Patients with Dukes Stage B and C Colorectal Cancer", \emph{Clin Cancer Res} \bold{15}(24):7642-7651.

Bild AH, Yao G, Chang JT, Wang Q, Potti A, Chasse D, Joshi MB, Harpole D, Lancaster JM, Berchuck A, Olson JA Jr, Marks JR, Dressman HK, West M, Nevins JR. (2006) "Oncogenic pathway signatures in human cancers as a guide to targeted therapies", \emph{Nature}, \bold{439}(7074):274-275.
}
\examples{
data(jorissen.colon.ras)
}
\keyword{ data }
