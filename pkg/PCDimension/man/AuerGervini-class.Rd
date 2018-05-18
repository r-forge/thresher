\name{AuerGervini-class}
\docType{class}
\alias{PCDimension}
\alias{AuerGervini-class}
\alias{plot,AuerGervini,missing-method}
\alias{summary,AuerGervini-method}
\alias{AuerGervini}
\alias{agDimension}
\title{Estimating Number of Principal Components Using the Auer-Gervini Method}
\description{
Auer and Gervini [1] described a graphical Bayesian method for
estimating the number of statistically significant principal
components.  We have implemented their method in the \code{AuerGervini}
class, and enhanced it by automating the final selection.
}
\usage{
AuerGervini(Lambda, dd=NULL, epsilon = 2e-16)
agDimension(object, agfun=agDimTwiceMean)
}
\arguments{
  \item{Lambda}{Either a \code{SamplePCA} object or a numerical vector
    of variances from a principal components analysis.}
  \item{dd}{A vector of length 2 containing the dimensions of the data
    used to created the Auer-Gervini object. If \code{Lambda} is a
    \code{SamplePCA} object, then the dimensions are taken from it,
    ignoring the \code{dd} argument.}
  \item{epsilon}{A numeric value.  Used to remove any variances that are
    less than \code{epsilon}; defaults to 2e-16.  Should only be needed in
    rare cases where negative variances show up because of round-off error.}
  \item{object}{An object of the \code{AuerGervini} class.}
  \item{agfun}{A function that takes one argument (a vector of step
    lengths) and returns a logical vector of the same length (where true
    indicates "long"  as opposed to "short" steps).}
}
\section{Objects from the Class}{
Objects should be created using the \code{AuerGervini} constructor.
}
\section{Slots}{
  \describe{
    \item{\code{Lambda}:}{A \code{numeric} vector containing
      the explained variances in decreasing order. }
    \item{\code{dimensions}}{Numeric vector of length 2 containing the
      dimnesions of the underlying data matrix.}
    \item{\code{dLevels}:}{Object of class \code{numeric}; see details }
    \item{\code{changePoints}:}{Object of class \code{numeric}; see details }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "AuerGervini", y = "missing")}: ... }
    \item{summary}{\code{signature(object = "AuerGervini")}: ... }
	 }
}
\details{
  The Auer-Gervini method for determining the number of principal
  components is based on a Bayesian model that assaerts that the
  vector of explained variances (eigenvalues) should have the form
  \deqn{a_1 \le a_2 \le \dots \le a_d < a_{d+1} = a_{d+2} = \dots a_n}
  with the goal being to find the true dimension \eqn{d}. They consider
  a set of prior distributions on \eqn{d \in \{1, \dots, n\}} that decay
  exponentially, with the rate of decay controlled by a parameter
  \eqn{\theta}.  For each value of \eqn{\theta}, one selects the value
  of \eqn{d} that has the maximum a posteriori (MAP) probability. Auer
  and Gervini show that the dimensions selected by this procedure write
  \eqn{d} as a non-increasing step function of \eqn{\theta}. The values
  of \eqn{\theta} where the steps change are stored in    the
  \code{changePoints} slot, and the corresponding \eqn{d}-values are 
  stored in the \code{dLevels} slot.

  Auer and Gervini go on to advise using their method as a graphical
  approach, manually (or visually?) selecting the highest step that is
  "long". Our implementation provides several different algorithms for
  automatically deciding what is "long" enough.  The simplest (but
  fairly naive) approach   is to take anything that is longer than twice
  the mean; other algorithms are described in
  \code{\link{agDimFunction}}.
}
\value{
  The \code{AuerGervini} function constructs and returns an object of
  the \code{AuerGervini} class.

  The \code{agDimension} function computes the number of significant
  principal components.  The general idea is that one starts by
  computing the length of each step in the Auer-Gerivni plot, and must
  then separate these into "long" and "short" classes.  We provide a
  variety of different algorithms to carry out this process; the
  default algorithm in the function \code{\link{agDimTwiceMean}} defines
  a step as "long" if it more than twice the mean step length.
}
\references{
  [1] P Auer, D Gervini.
  Choosing principal components: a new graphical method based on Bayesian model selection.
  Communications in Statistics-Simulation and Computation 37 (5),
  962-977.

  [2] Wang M, Kornbla SM, Coombes KR.
  Decomposing the Apoptosis Pathway Into Biologically Interpretable
  Principal Components.
  Preprint: bioRxiv, 2017. <doi://10.1101/237883>.

}
\author{Kevin R. Coombes <krc@silicovore.com>}
\seealso{
  \code{\link{agDimFunction}} to get a complete list of the functions
  implementing different algorithms to separate the step lengths into
  two classes.
}
\examples{
showClass("AuerGervini")
# simulate variances
lambda <- rev(sort(diff(sort(c(0, 1, runif(9))))))
# apply the Auer-Gervini method
ag <- AuerGervini(lambda, dd=c(3,10))
# Review the results
summary(ag)
agDimension(ag)
agDimension(ag, agDimKmeans)
# Look at the results graphically
plot(ag, agfun=list(agDimTwiceMean, agDimKmeans))
}
\keyword{classes}
\keyword{models}
