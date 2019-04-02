\name{removeDuplicateFeatures}
\alias{removeDuplicateFeatures}
\title{
Remove Duplicate Features from a Binary Matrix Object
}
\description{
The \code{removeDuplicateFeatures} function removes duplicate columns from 
a \code{binaryMatrix} object in the \code{Mercator} package.
}
\usage{
removeDuplicateFeatures(binmat)
}
\arguments{
  \item{binmat}{An object of class \code{binaryMatrix}.}
}
\details{
In some analyses, it may be desirable to remove duplicate features to collapse 
a group of identical, related events to a single feature, to prevent overweighting 
when clustering.

Removal of duplicate features is not required for performance of the 
\code{binaryMatrix} or \code{Mercator} objects and associated functions.

The \code{@history} of the \code{binaryMatrix} object documents removal of 
duplicate features.

Future versions of this package may include functionality to store the identities 
of duplicate features removed.
}
\value{
Returns an object of class \code{binaryMatrix} with duplicate columns removed.
}
\author{Kevin R. Coombes <krc@silicovore.com}
\note{
Transposing the \code{binaryMatrix} can allow the \code{removeDuplicateFeatures} 
function to be applied to both features and observations, if desired.

Features containing exclusively 0s or 1s may interfere with performance of 
\code{removeDuplicateFeatures}.
}
\examples{my.binmat <- removeDuplicateFeatures(my.binmat)}
\keyword{ cluster }
