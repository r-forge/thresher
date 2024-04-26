\name{Thresher-data}
\alias{Thresher-data}
\alias{sigma}
\alias{savedSims}
\alias{savedReap}
\docType{data}
\title{Thresher and Reaper Simulated Data}
\description{
  This data set contains five related simulated data sets, along with
  the Thresher and Reaper objects obtained by clustering the data sets
  and removing outliers. 
}
\usage{
data(savedSims)
}
\format{
  \describe{
  \item{\code{sigma}}{: A list of length five; each entry is a covariance matrix
    used to simulate data.}
  \item{\code{savedSims}}{: A list of length five; each entry is a
    \code{\link{SimThresher}} object obtained by simulating data from
    one of the covariance matrices and running the \code{\link{Thresher}}
    algorithm.}
  \item{\code{savedReap}}{: A list of length five; each entry is a
    \code{\link{Reaper}} object obtained by applying the \code{Reaper}
    function.}
  }
} 
\source{
  The simulated data sets were prepared by running the script
  \code{Examples/makeSims.R} that is installed along with the
  \code{Thresher} package. The five covariance matrices vary in the
  number of correlated subgroups (one or two) and whether they include
  both positively and negatively correlated variables, or just
  positively correlated ones. Each data set also includes two
  uncorrelated "noise" variables that should automatically be removed by
  the Reper-Thresher algorithms.
}
\references{
  Wang M, Abrams ZB, Kornblau SM, Coombes KR.
  Thresher: determining the number of clusters while removing outliers.
  BMC Bioinformatics, 2018; 19(1):1-9. doi://10.1186/s12859-017-1998-9.
}
\keyword{datasets}

