\name{eucdist}
\alias{eucdist}
\docType{data}
\title{An example Euclidean distance matrix}
\description{
  The Euclidean distance matrix between 300 objects, used to illustrate
  the SillyPutty algorithms.
}
\usage{data(eucdist)}
\format{
  A \code{dist} object representing Euclidean distances.
}
\details{
  This data set was generated in the \code{SillyPutty} vignette from the
  tools in the \code{Umpire} R package.
}
\source{
  This data set was generated in the \code{SillyPutty} vignette from the
  tools in the \code{Umpire} R package, and saved using code that is no
  hidden and disabled in the vignette source.
}
\examples{
data(eucdist)
class(eucdist)
attr(eucdist, "Size")
}
\keyword{datasets}
