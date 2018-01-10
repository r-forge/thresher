\name{unitize}
\alias{unitize}
\title{
 Convert a Vector to Unit Length
}
\description{
Rescales each column of a matrix to produce vectors of length one.
}
\usage{
unitize(mat)
}
\arguments{
  \item{mat}{A matrix of real numbers.}
}
\details{
No details beyond the simple description are requires; it is implemented
exactly the way you would suspect.
}
\value{
  A matrix of the same size as the input matrix.
}
\author{
Kevin R. Coombes <krc@silicovore.com>
}
\examples{
simmat <- matrix(rnorm(12), 3, 4)
U <- unitize(simmat)
apply(U^2, 2, sum)  # check unit length
simmat/U            # view normalization factors
}
\keyword{array }
