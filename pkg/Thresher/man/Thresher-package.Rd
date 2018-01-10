\name{Thresher-package}
\alias{Thresher-package}
\docType{package}
\title{
Threshing and Reaping for Principal Components
}
\description{
  The \code{Thresher} package defines classes used to identify outliers
  (threshing) and to compute the number of significant principal
  components and number of clusters (reaping) in a joint application of
  PCA and model-based clustering.
}
\details{
\tabular{ll}{
Package: \tab Thresher\cr
Type: \tab Package\cr
Version: \tab 0.11.5\cr
Date: \tab 2018-01-09\cr
License: \tab Apache 2.0\cr
Depends: \tab ClassDiscovery, PCDimension\cr
Imports: \tab methods, stats, graphics, MASS, colorspace, movMF, ade4, oompaBase\cr
}
Standard usage is to apply the \code{Thresher} function to a data set in
order to estimate the principal component dimension and identifty
outliers. You then apply the \code{Reaper} function to actually remove the
outliers and assign the remaining objects to clusters.
}
\author{
Kevin R. Coombes and Min Wang

Maintainer:  <krc@silicovore.com>
}
\references{
  Wang M, Abrams ZA, Kornblau SM, Coombes KR.
  Thresher: detemining the number of clusters while removing outliers.
  BMC Bioinformatics. 2018. In press.
}
\seealso{
  \code{\link{PCDimension}}
}
\keyword{ package }
