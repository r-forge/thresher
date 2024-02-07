\name{HCSP}
\alias{HCSP}
\title{Combining Hierarchical Clustering with SillyPutty}
\description{
  Our simulations revealed that the fastest and most accuirate
  clustering algorithm for modest-sized contiuous data sets is the
  combination of hierarchical clustering (with Ward's linkage rule)
  followed by SillyPutty. The function \code{HCSP} implements this
  combination.
}
\usage{
  HCSP(dis, K, method = "ward.D2", \dots)
}
\arguments{
  \item{dis}{An object of class \code{dist} representing a distance matrix.}
  \item{K}{The desired number of clusters.}
  \item{method}{Sane as the corresponding argument for \code{\link{hclust}}.
    We recommend not changing it.}
  \item{\dots}{Extra arguments to the \code{SillyPutty} function.}
}
\value{
  A list containing two items: \code{hc}, the results of hierarchical
  clustering, and \code{sp}, a \code{SillyPutty} object by applying the
  algorithm to the result of cutting the dendrogram to produce \code{K}
  clusters. 
}
\details{
  The \code{HCSP} function that first runs hierarchical clustering, then
  applies the \code{SillyPutty} algorithm.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com}
}
\references{
  Polina Bombina, Dwayne Tally, Zachary B. Abrams, Kevin R. Coombes.
  SillyPutty: Improved clustering by optimizing the silhouette width,
  bioRxiv 2023.11.07.566055; doi: https://doi.org/10.1101/2023.11.07.566055
}
\examples{
data(eucdist)
set.seed(1234)
twostep <- HCSP(eucdist, K=5)
sw <- cluster::silhouette(twostep$sp@cluster, eucdist)
plot(sw)
}
\keyword{
  cluster
}
