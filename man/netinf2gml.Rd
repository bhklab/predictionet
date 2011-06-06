\name{netinf2gml}
\alias{netinf2gml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to create an \code{\link[igraph]{igraph}} object and export a network to a GML readable by Cytoscape
}
\description{
This function creates, from a network inferred from \link[predictionet]{netinf} or \link[predictionet]{netinf.cv}, an \code{\link[igraph]{igraph}} object and export this network to a GML readable by Cytoscape.
}
\usage{
netinf2gml(object, edge.info, node.info, file = "predictionet")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object returns by \code{netinf} or \code{netinf.cv}
}
  \item{edge.info}{
matrix of values representing the statistics for each edge; parents in rows, children in columns. A list of matrices could be provided, names of the list will then be used to describe the statistics in Cytoscape
}
	\item{node.info}{
vector of values representing the statistics for each node; parents in rows, children in columns. A list of vectors could be provided, names of the list will then be used to describe the statistics in Cytoscape
}
	\item{file}{
name of the GML file to be saved.
}
}
\details{
The GML file created by this function has been tested on Cytoscape 2.8.1; a Vizmap property file of the same name is also created and could be imported into Cytoscape ("preditionet_vizmap2") so the information for each node and edge are displayed correctly.
}
\value{
an \link[igraph]{igraph} object
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
\\code{RCytoscape}
}
\examples{
## load gene expression data for colon cancer data, list of genes related to RAS signaling pathway and the corresponding priors
data(expO.colon.ras)
## number of genes to select for the analysis
genen <- 10
## select only the top genes
goi <- dimnames(annot.ras)[[1]][order(abs(log2(annot.ras[ ,"fold.change"])), decreasing=TRUE)[1:genen]]
mydata <- data.ras[ , goi, drop=FALSE]
myannot <- annot.ras[goi, , drop=FALSE]
mypriors <- priors.ras[goi, goi, drop=FALSE]
mydemo <- demo.ras
## infer global network from data and priors
mynet <- netinf.cv(data=mydata, categories=3, priors=mypriors, priors.count=TRUE, priors.weight=0.5, maxparents=3, method="regrnet", nfold=3, seed=54321)

## create an igraph obkect and export it into a GML file
netinf2gml(object=mynet, file = "predictionet")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graph }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
