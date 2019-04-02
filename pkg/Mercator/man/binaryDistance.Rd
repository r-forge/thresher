\name{binaryDistance}
\alias{binaryDistance}
\alias{jaccardSimilarity}
\alias{jaccardDistance}
\alias{sokalMichenerSimilarity}
\alias{sokalMichenerDistance}
\alias{hammingDistance}
\alias{russellRaoSimilarity}
\alias{russelRaoDistance}
\alias{pearsonSimilarity}
\alias{pearsonDistance}
\alias{goodmanKruskalSimilarity}
\alias{goodmanKruskalDistance}
\title{
Binary Distance Definitions Across Methods
}
\description{
The \code{binaryDistance} function defines a various similarity or distance
measures between binary vectors, which represent the first step in the 
algorithm underlying the /code{Mercator} visualizations.
}
\usage{
binaryDistance(mat, metric)
}
\arguments{
  \item{mat}{An object of class \code{matrix}.}
  \item{metric}{An object of class \code{character} limited to the names of 
  10 selected distance metrics: \code{jaccard}, \code{sokalMichener}, \code{hamming},
  \code{russellRao}, \code{pearson}, \code{goodmanKruskal}, \code{manhattan},
  \code{canberra}, \code{binary}, or \code{euclid}.
  }
}
\value{

}
\details{
Similarity or difference between binary vectors can be calculated using a
variety of distance measures. In the main reference (below), Choi and colleagues
reviewed 76 different measures of similarity of distance between binary vectors.
They also produced a hierarchical clustering of these measures, based on the 
correlation between their distance values on multiple simulated data sets. For
metrics that are highly similar, we chose a single representative.

Cluster 1, represented by the \code{jaccard} distance, contains Dice & Sorenson, Ochiai, 
Kulcyznski, Bray & Curtis, Baroni-Urbani & Buser, and Jaccard.

Cluster 2, represented by the \code{sokalMichener} distance, contains Sokal & Sneath, 
Gilbert & Wells, Gower & Legendre, Pearson & Heron, Hamming, and Sokal & Michener. 
Also within this cluster are 4 distances represented independently within this function: 
\code{hamming}, \code{manhattan}, \code{canberra}, and \code{euclidean} distances

Cluster 3, represented by the \code{russellRao} distance, contains Driver & Kroeber,
Forbes, Fossum, and Russell & Rao.

The remaining metrics are more isolated, without strong clustering. We considered a few 
examples, including the Pearson distance (\code{pearson}) and the Goodman & Kruskal distance 
(\code{goodmanKruskal}). The \code{binary} distance is also included.
}
\references{
Choi SS, Cha SH, Tappert CC, A Survey of Binary Similarity and Distance Measures. 
Systemics, Cybernetics, and Informatics. 2010; 8(1):43-48.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\note{The \code{labelAccuracy} function should probably allow the user
  to supply a list of distance metrics instead of relying on the
  hard-coded list internally.}
\seealso{
  Hierarchical clustering is implemented in the \code{\link{hclust}}
  function.  We use the extended set of distance metrics provided by the
  \code{\link{distanceMatrix}} function from the ClassDiscovery package.
  This set includes all of the  metrics from the \code{\link{dist}}
  funciton. 
}
\examples{
factor1 <- sample(c("A", "B", "C"), 30, replace=TRUE)
factor2 <- rep(c("X", "Y", "Z"), each=10)
tab <- table(factor1, factor2)
matchLabels(tab)
labelMatcher(tab)
R <- remap(factor1, factor2)
table(R, factor2) # remapping
table(R, factor1) # cross-comparison
}
\keyword{ cluster }
\keyword{ classif }
\keyword{ multivariate }