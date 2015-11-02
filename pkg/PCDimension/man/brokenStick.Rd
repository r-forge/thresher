\name{brokenStick}
\alias{brokenStick}
\alias{bsDimension}
\title{
  The Broken Stick Method
}
\description{
  The Broken Stick model is one proposed method for estimating the
  number of statistically significant principal components. 
}
\usage{
brokenStick(k, n)
bsDimension(lambda, FUZZ = 0.005)
}
\arguments{
  \item{k}{An integer between 1 and \code{n}.}
  \item{n}{An integer; the total number of principal components.}
  \item{lambda}{The set of variances from each component from a principal
    components analysis.  These are assumed to be already sorted in
    decreasing order. You can also supply a \code{\link{SamplePCA}}
    object, and the variances will be automatically extracted.}
  \item{FUZZ}{A real number; anything smaller than \code{FUZZ} is
    assumed to equal zero for all practical purposes.}
}
\details{
  The Broken Stick model is one proposed method for estimating the
  number of statistically significant principal components.  The idea is
  to model \eqn{N} variances by taking a stick of unit length and breaking it
  into \eqn{N} pieces by randomly (and simultaneously) selecting break
  points from a uniform distribution.
}
\value{
  The \code{brokenStick} function returns, as a real number, the
  expected value of the \code{k}-th longest piece when breaking a
  stick of length one into \code{n} total pieces.  Most commonly used
  via the idiom \code{brokenStick(1:N, N)} to get the entire vector of
  lengths at one time.

  The \code{bsDimension} function returns an integer, the number of
  significant components under this model.  This is computed by finding
  the last point at which the observed variance is bugger than the
  expected value under the broken stick model by at least \code{FUZZ}.
}
\references{
Jackson, D. A. (1993). Stopping rules in principal components analysis:
a comparison of heuristical and statistical approaches. Ecology 74,
2204--2214. 

Legendre, P. and Legendre, L. (1998) Numerical Ecology. 2nd English ed. Elsevier.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\seealso{
  Better methods to address this question are based on the Auer-Gervini
  method; see \code{\link{AuerGervini}}.
}
\examples{
brokenStick(1:10, 10)
sum( brokenStick(1:10, 10) )
fakeVar <- c(30, 20, 8, 4, 3, 2, 1)
bsDimension(fakeVar)
}
\keyword{ models }
\keyword{ cluster }
