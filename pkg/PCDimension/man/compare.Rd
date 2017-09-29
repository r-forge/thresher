\name{compareAgDimMethods}
\alias{compareAgDimMethods}
\title{
  Compare Methods to Divide Steps into "Long" and "Short"
}
\description{
  Auer and Gervini developed a Bayesian graphical method to determine
  the number \eqn{d} of significant principal components; a brief
  overview is included in the help for the \code{\link{AuerGervini}}
  class.  The output of their method is a step function that displays
  the maximum a posteriori (MAP) choice of \eqn{d} as a step function of
  a one-parameter family of prior distributions, and they recommend
  choosing the highest "long" step.  The functions described here help
  automate the process of dividing the step lengths into "long" and
  "short" classes.
}
\usage{
  compareAgDimMethods(object, agfuns)
}
\arguments{
  \item{object}{An object of the \code{\link{AuerGervini}} class}
  \item{agfuns}{A list of functions}
}
\details{
  This method simply iterates over the list of functions that implement
  different algorithms/methods to determine the PC dimension.
}
\value{
  Returns an integer vector of te same length as the list of
  \code{agfuns}, containing the number of significant principal
  components computed by each method.
}
\references{
  P Auer, D Gervini.
  Choosing principal components: a new graphical method based on Bayesian model selection.
  Communications in Statistics-Simulation and Computation 37 (5), 962-977
}
\author{
  Kevin R. Coombes <krc@silicovore.com>, Min Wang <wang.1807@osu.edu>.
}
\seealso{
  \code{\link{AuerGervini}}, \code{\link{agDimension}}.
}
\examples{
# simulate variances
lambda <- rev(sort(diff(sort(c(0, 1, runif(9))))))
# apply the Auer-Gervini method
ag <- AuerGervini(lambda, dd=c(3,10))
# try different methods
agfuns <- list(twice=agDimTwiceMean,
               km=agDimKmeans,
               cpt=agDimCPT)
compareAgDimMethods(ag, agfuns)
}
\keyword{models}
