\name{eucdist}
\alias{eucdist}
\alias{trueGroups}
\docType{data}
\title{An example Euclidean distance matrix}
\description{
  The Euclidean distance matrix between 300150 objects, used to illustrate
  the SillyPutty algorithms.
}
\usage{data(eucdist)}
\format{
  The binary R data file ocntains two objects, First, a \code{dist}
  object representing Euclidean distances between 150 samples. Second, a
  vector of the known (simualted) true groups to which each sample
  belongs.
}
\details{
  This data set was generated in the \code{SillyPutty} vignette from
  tools in the \code{Umpire} R package. The simulated data was intended
  to have five different clusters, all of approximately the same size.
  Noise ws added to make the clusters somehwat harder for most
  algorithms to distinguish.TRhe same data set is used in most of the
  examples in the man pages.
}
\source{
  This data set was generated in the \code{SillyPutty} vignette from the
  tools in the \code{Umpire} R package, and saved using code that is now
  hidden and disabled in the vignette source.
}
\examples{
data(eucdist)
class(eucdist)
attr(eucdist, "Size")
}
\keyword{datasets}
