\name{CytobandData}
\alias{RcytoGPS}
\title{CytobandData}
\description{
  TThe purpose of this function is to convert a frequency dataframe with the CL columns (chromosome, loc.start, loc.end, the band, and the stain) into a s4 object that can be used with RCytoGPS function images to generate a multitude of graphs to give you a clearer idea of what your data is expressing.
}
\usage{
CytobandData(data)
}
\arguments{
  \describe{
    \item{\code{data}}{ A dataframe that has numeric values as well as CL columns.}
  }
}
\value{
 Returns an s4 object that can be used to generate RCytoGPS images.

}
\source{

}
\references{

}
\examples{
temp <- readLGF(folder = wd)
cytoData <- data.frame(temp[["CL"]],temp[["frequency"]])
bandData <- CytobandData(cytoData)

}
\keyword{datasets}
