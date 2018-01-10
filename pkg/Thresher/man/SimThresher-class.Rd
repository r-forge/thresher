\name{SimThresher-class}
\Rdversion{1.1}
\docType{class}
\alias{SimThresher-class}
\alias{SimThresher}
\alias{image,SimThresher-method}
\alias{makeFigures,SimThresher-method}

\title{Class \code{"SimThresher"}}
\description{
The \code{SimThresher} class is used to simulate \code{\link{Thresher}} objects.
}
\section{Objects from the Class}{
Objects should be created using the \code{SimThresher} constructor.
}
\usage{
SimThresher(ss, nSample, nm = deparse(substitute(ss)), rho = NULL,
            agfun = agDimTwiceMean, ...)
}
\arguments{
  \item{ss}{A covariance matrix.}
  \item{nSample}{An integer; the number of samples to simulate.}
  \item{nm}{A character string; the name of this object.}
  \item{rho}{A numeric vector; the correlation between different
    variables. If NULL, then these are obtained from the covariance
    matrix.}
  \item{agfun}{A function used by the \code{\link{AuerGervini}}
    function to determine the number of principal components.}
  \item{\dots}{Parameters to be passed to the \code{\link{Thresher}}
    constructor.}
}
\section{Slots}{
  \describe{
    \item{\code{nSample}:}{An integer; the number of simulated samples.}
    \item{\code{covariance}:}{A covariance matrix.}
    \item{\code{rho}:}{A vector of correlation coefficients; essentially
      the unique values in the upper triangular part of the covariance
      matrix.}
  }
}
\section{Extends}{
Class \code{\linkS4class{Thresher}}, directly.
}
\section{Methods}{
  \describe{
    \item{image}{\code{signature(x = "SimThresher")}: Produces an image
      of the covariance matrix.}
    \item{makeFigures}{\code{signature(object = "SimThresher")}: This is a
      convenience function to produce a standard set of figures. In
      addition tot he plots preodcued for\code{Thresher} object, this
      function also produces an image of te covariance matrix used in
      the simulations.
      If the \code{DIR} argument is
      non-null, it is treated as the name of an existing directory where the
      figures are stored as PNG files. Otherwise, the figures are
      displayed interactively, one at a time, in a window on screen.}
  }
}
\details{
  Basically, given a number of samples and a covariance matrix,
  simulate a data matrix of the appropriate size and multivariate normal
  structure by assuming that all of the means are zero. After simulating
  the data, we apply the Thresher algorithm. The result is an
  object that combines the simulation parameters, simulated data, and
  fitted model.
}
\value{
  The \code{SimThresher} function returns an object of the
  \code{SimThresher} class.
}
\references{
  Wang M, Abrams ZB, Kornblau SM, Coombes KR.
  Thresher: determining the number of clusters while removing outliers.
  BMC Bioinformatics, 2018; 19(1):1-9. doi://10.1186/s12859-017-1998-9.
}
\author{
 Kevin R. Coombes <krc@silicovore.com>, Min Wang.
}
\seealso{
  \code{\link{Thresher}}, \code{\link{Reaper}}.
}
\examples{
set.seed(250264)
rho <- 0.5
nProtein <- 16
splinter <- sample((nProtein/2) + (-3:3), 1)
sigma1 <- matrix(rho, ncol=nProtein, nrow=nProtein)
diag(sigma1) <- 1
st <- SimThresher(sigma1, nSample=300)
image(st, col=redgreen(64), zlim=c(-1,1))
screeplot(st, col='pink', lcol='red')
plot(st)
scatter(st)
heat(st)
}
\keyword{classes}
