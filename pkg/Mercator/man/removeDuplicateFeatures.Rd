\name{removeDuplicates}
\alias{removeDuplicates}
\alias{removeDuplicateFeatures}
\title{
Remove Duplicate Features or Samples from a Binary Matrix Object
}
\description{
The \code{removeDuplicates} function removes duplicate columns from 
a \code{binaryMatrix} object in the \code{Mercator} package.
}
\usage{
removeDuplicates(object)
removeDuplicateFeatures(object)
}
\arguments{
  \item{object}{An object of class \code{binaryMatrix}.}
}
\details{
In some analyses, it may be desirable to remove duplicate features to collapse 
a group of identical, related events to a single feature, to prevent overweighting 
when clustering. Historically, this funciton was called
\code{removeDuplicateFeatures}. That name is still retained for
backwards compatibility, but it may be deprecated in future versions in
favor of \code{removeDuplicates}.
In the same way, for some clustering applications,
it may be usedful to remove duplicate samples, or those that have an
identical feature set. 

Removal of duplicates is not required for performance of the 
\code{binaryMatrix} or \code{Mercator} objects and associated functions.

The \code{history} slot of the \code{binaryMatrix} object documents the removal of 
duplicates.

Future versions of this package may include functionality to store the identities 
of any duplicates that were removed.
}
\value{
Returns an object of class \code{binaryMatrix} with duplicate columns removed.
}
\author{Kevin R. Coombes <krc@silicovore.com>, Caitlin E. Coombes}
\note{
Transposing the \code{binaryMatrix} can allow the \code{removeDuplicates} 
function to be applied to both features and observations, if desired.

Features containing exclusively 0s or 1s may interfere with performance of 
\code{removeDuplicates}.
}
\examples{
my.matrix <- matrix(rbinom(50*100, 1, 0.15), ncol=50)
my.matrix <- cbind(my.matrix, my.matrix[, 1:5]) # add duplicates
dimnames(my.matrix) <- list(paste("R", 1:100, sep=''),
                            paste("C", 1:55, sep=''))
my.binmat <- BinaryMatrix(my.matrix)
dim(my.binmat)
my.binmat <- removeDuplicates(my.binmat)
dim(my.binmat)
}
\keyword{ cluster }
