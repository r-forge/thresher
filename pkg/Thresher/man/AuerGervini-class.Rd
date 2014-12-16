\name{AuerGervini-class}
\Rdversion{1.1}
\docType{class}
\alias{AuerGervini-class}
\alias{plot,AuerGervini,missing-method}
\alias{summary,AuerGervini-method}
\alias{AuerGervini}
\alias{agDimension}
\title{Estimating Number of Principal Components Using the Auer-Gervini Method}
\description{
Auer and Gervini [1] described a graphical Bayesian method for
estimating the number of statisticailly significant principal
components.  We have implemented their method in the \code{AuerGervini}
class, and enhanced it by automating the final selection.
}
\usage{
AuerGervini(Lambda, epsilon = 2e-16)
agDimension(object)
}
\arguments{
  \item{Lambda}{Either a \code{SamplePCA} object or a numerical vector
    of variances from a principal components analysis.}
  \item{epsilon}{A numeric value.  Used to remove any variances that are
  less than \code{epsilon}; defaults to 2e-16.  Should only be needed in
  rare cases where negative variances show up because of round-off error.}
  \item{object}{An object of the \code{AuerGervini} class.}
}
\section{Objects from the Class}{
Objects should be created using the \code{AuerGervini} constructor.
}
\section{Slots}{
  \describe{
    \item{\code{Lambda}:}{Object of class \code{"numeric"}, containing
      the vector of explained variances. }
    \item{\code{dLevels}:}{Object of class \code{"numeric"}; see details }
    \item{\code{changePoints}:}{Object of class \code{"numeric"}; see details }
    \item{\code{lowerBounds}:}{Object of class \code{"numeric"}; see details  }
    \item{\code{upperBounds}:}{Object of class \code{"numeric"}; see details  }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "AuerGervini", y = "missing")}: ... }
    \item{summary}{\code{signature(object = "AuerGervini")}: ... }
	 }
}
\details{
 TODO: Explain the Auer-Gervinii computationsd, and how they contribute
 to the slots of this class.	 
}
\value{
  The \code{AuerGervini} function constructs and returns an object of
  the \code{AuerGervini} class.

  The \code{agDimension} function computes the number of significant
  principal components.
}
\references{
  [1] P Auer, D Gervini.
  Choosing principal components: a new graphical method based on Bayesian model selection.
  Communications in Statistics-Simulation and Computation 37 (5), 962-977
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\seealso{
}
\examples{
showClass("AuerGervini")
# simulate variances
lambda <- rev(sort(diff(sort(c(0, 1, runif(9))))))
# apply the Auer-Gervini method
ag <- AuerGervini(lambda)
# Review the results
summary(ag)
agDimension(ag)
# Look at the results graphically
plot(ag)
}
\keyword{classes}
