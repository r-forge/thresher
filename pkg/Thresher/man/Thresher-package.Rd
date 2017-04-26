\name{Thresher-package}
\alias{Thresher-package}
\alias{Thresher}
\docType{package}
\title{
What the package does (short line)
Threshing and Reaping for Principal Components
}
\description{
  The Thresher package defines the classes used to identify
  outliers (threshing) and compute the number of significant principal
  components and number of clusters (reaping) in a joint application
  of PCA and hierarchical clustering.
}
\details{
\tabular{ll}{
Package: \tab Thresher\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-11-20\cr
License: \tab What license is it under?\cr
Depends: \tab methods\cr
}
Standard usage is to apply the \code{Thresher} function to a data set in
order to estimate the princiapl component dimension and identifty
outliers. You then apply the \code{Reaper} function to actually remove the
outliers and asssign the remaining objects to clusters.
}
\author{
Kevin R. Coombes and Min Wang

Maintainer:  <krc@silicovore.com>
}
\references{
  Submitted to some jounral.
}
\keyword{ package }
\seealso{
  \code{\link[PCDimension:PCDimesnion-package]{PCDimension}}
}
\examples{
}
