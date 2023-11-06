\name{RandomSillyPutty-class}
\docType{class}
\alias{RandomSillyPutty-class}
\alias{RandomSillyPutty}
\alias{plot,RandomSillyPutty,matrix-method}
\alias{plot,RandomSillyPutty,missing-method}
\alias{summary,RandomSillyPutty-method}
\alias{plot,RSPSummary,missing-method}
\alias{plot,RSPSummary,matrix-method}


\title{Running SillyPutty From Multiple Random Initial Clusterings}
\description{
  A function to perform cluster assignments on a distance matrix based
  on optimizing silhouette width. The cluster assignments are based on
  maximum and minimum silhouette width scores obtained from N iterations.
}
\usage{
RandomSillyPutty(distobj, K, N = 100, verbose = FALSE, \dots)
\S4method{plot}{RandomSillyPutty,matrix}(x, y, distobj, col = NULL, \dots)
\S4method{plot}{RandomSillyPutty,missing}(x, y, \dots)
\S4method{summary}{RandomSillyPutty}(object, \dots)
\S4method{plot}{RSPSummary,missing}(x, y, \dots)
}
\arguments{
  \item{distobj}{An object of class \code{\link{dist}}.}
  \item{K}{The number of clusters.}
  \item{N}{The number of iterations you want to run.}
  \item{verbose}{A logical value; should you print info while working}
  \item{\dots}{Extra arguments to the \code{SillyPutty} function or to
    generic methods.}
  \item{x}{An object of the \code{RandomSillyPutty} or \code{RSPSummary}
    class.}
  \item{object}{An object of the \code{RandomSillyPutty} class.}
  \item{y}{A layout matrix.}
  \item{col}{A character vector containing color names.}
}

\section{Slots}{
  \describe{
    \item{\code{MX}:}{An integer vector containing cluster assignment that had
      the best silhouette width from running the iterations}
    \item{\code{MN}:}{An integer vector containing cluster assignment that had
      the worst silhouette width from running the iterations}
    \item{\code{ave}:}{An integer vector of average distribution of the 
      silhouette width scores throughout N iterations}
    \item{\code{labels}:}{A data frame of the cluster assignments of the best 
      silhouette width score.}
    \item{\code{minItemSW}:}{A list of silhouette width scores of all N 
      iterations}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "RandomSillyPutty", y = "matrix")}: Plot the
  clusterings with the maximum and minimum global silhouette widths.}
    \item{summary}{\code{signature(x = "RandomSillyPutty")}: .}
  }
}
\value{
  The constructor fnuction, \code{RandomSillyPutty}, returns an object
  of the \code{RandomSillyPutty} class.
}
\details{
  The \code{RandomSillyPutty} function reads and processes one distance matrix
  at a time, with a precomputed cluster number, and a number N iterations. 
  RandomSillyPutty returns an s4 object. The \code{MX} component of this 
  structure contains an integer vector that has a cluster assignment based on the
  best scoring silhouette width score from N iterations. The \code{MN} contains 
  an integer vector that has a cluster assignment based on the worst scoring 
  silhouette score from N iterations. The \code{ave} contains the average 
  silhouette width of all N iteration. The \code{labels} is a dataframe containing
  the cluster assignment of the best scoring slihouette width score per iteration.
  The \code{minItemSW} is a list containing the silhouette width score of all N 
  iterations.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  Dwayne G. Tally \email{dtally110@hotmail.com}
}

\references{
  Pending.
}
\examples{
data(eucdist)
# 'eucdist' is the Euclidean distane matrix (i.e., 'dist' object) from
# a simulated data set with 500 elements and 5 groups.
set.seed(12)
y <- RandomSillyPutty(eucdist, 6, N = 100)
summary(y)
}
\keyword{
  cluster
}
