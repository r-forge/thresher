\name{findClusterNumber}
\alias{findClusterNumber}
\title{Using SillyPutty to find the number of clusters}
\description{
  A function that is designed to find an approximation of the true
  number. K, of clusters in a dataset. the \code{findClusterNumber}
  function calls \code{RandomSillyPutty} for each value of K in the
  range fro \code{start} to \code{end}, performing \code{N} random
  starts eacyh time.

  NOTE: start must be > 1, and the function can be slow depending on how
  complex the dataset is and the number of N iterations. 
}
\usage{
  findClusterNumber(distobj, start,end, N = 100)
}
\arguments{
  \item{distobj}{An object of class {distance matrix}}
  \item{start}{The minimum cluster number for the range of clusters}
  \item{end}{The maximum cluster number for the range of clusters}
  \item{N}{Number of iterations}
}

\value{
  List:
    \itemize{
      \item{\code{list:} {A list containing the maximum silhouette width values per K clusters.}
      }
    }
  }
\details{
  The \code{findClusterNumber} function processes one distance matrix at
  a time, through N iterations. It returns a list. The \code{list} is a
  list of the maximum silhoutte width values obtained from N iterations 
  with their associated cluster number.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  Dwayne G. Tally \email{dtally110@hotmail.com}
}
\references{
  
}
\examples{
data(spExample)
set.seed(12)
y <- findClusterNumber(eucdist, start = 3, end = 7)
plot(names(y), y, xlab = "K", ylab = "Mean Silhouette Width",
     type = "b", lwd = 2, pch = 16)
}
\keyword{
  cluster
}
