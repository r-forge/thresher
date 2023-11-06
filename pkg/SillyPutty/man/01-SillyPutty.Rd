\name{SillyPutty-class}
\docType{class}
\alias{SillyPutty-class}
\alias{SillyPutty}
\title{Running SillyPutty}
\description{
  A function that takes in an already existing starting location based on
  unsupervised clustering attempts. I.G. Kmeans or Hieriarchical cluster
  assignment. SillyPutty optimizes the pre-exisitng cluster assignments based on
  the best silhouette width score.
}
\usage{
  SillyPutty(labels, dissim, maxIter = 1000, loopSize = 15, verbose = FALSE)
}
\arguments{
  \item{labels}{A numeric vector containing pre-computed cluster labels}
  \item{dissim}{An object of class \code{dist}; that is, a distance matrix.}
  \item{maxIter}{A numneric vetor of length one; the maximum number of
    individual steps, each of which reclassifies only one object}
  \item{loopSize}{How many steps to retain in momry to test if you have
    entered an infinite loop while rearranging objects.}
  \item{verbose}{A logical vector of length one; should you output a lot
    of information while running?}
}
\section{Slots}{
  \describe{
    \item{\code{cluster}:}{A list containing the adjusted cluster
      assignment that had the best silhouette width.}
    \item{\code{silhouette}:}{A dataframe containing the silhouette
      width scores.}
    \item{\code{minSW}:}{A silhouette double vector that contains the
      positive and negative version of the minimum silhouette width value.}
    \item{\code{meanSW}:}{A double vector that contains the average
      silhouette width value.}
  }
}
\value{The constructor function \code{SillyPutty}, returns an object of
  the \code{SillyPutty} class.
}
\details{
  The \code{SillyPutty} function processes a pre-computed cluster assignment
  along with a distance metric and returns a s4 object. The \code{cluster} 
  component is a list of the new cluster assignments based on best 
  silhouette width score. The \code{silhouette} is a dataframe containing the 
  silhouette width score calculated by SillyPutty. The \code{minSw} contains 
  a positive and negative version of the minimum silhouette width score.
  The \code{meanSW} is a double vector that shows the average silhouette width
  score after applying SillyPutty to the cluster assignment.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  Dwayne G. Tally \email{dtally110@hotmail.com}
}
\references{
  Pending
}
\examples{
data(eucdist)
set.seed(12)
hc  <- hclust(eucdist, "ward.D2")
clues <- cutree(hc, k = 5)
hcSilly <- SillyPutty(clues, eucdist)
}
\keyword{
  cluster
}
  
