
\name{Idioformat}
\alias{RcytoGPS}
\title{Idioformat}
\description{
The purpose of this function is to generate and format the frequency based off of the LGF dataframe output from extractOneLGF. Idioformat takes the dataset from the LGF dataframe and column means the dataset, based off what cytogenetic event that is taking place. CL is the cytobandLocation object which is found within the package RCytoGPS that gives additional columns needed for frequency to generate images, and Idioformat also produces the size of the dataset so you can see the expected size of the files that are being read to aid in determining if readLGF is taking in the right input files.
}
\usage{
Idioformat(df, CL)
}
\arguments{
  \describe{
    \item{\code{df}}{A dataframe that has the numeric columns of the LGF model.}
    \item{\code{CL}}{A dataframe that is obtained within the RCytoGPS package. It contains the columns for chromosome, loc.start, loc.end, the band, and the stain.}
  }
}
\value{
   Returns a dataframe of the binary LGF model column mean, as well as the size of the dataset.
}
\source{

}
\references{

}
\examples{
temp <- readLGF(folder = wd)
 test <- lapply(temp[["raw"]], function(R) {
  Idioformat(R$LGF, CL = cytobandLocations)
 })

}
\keyword{datasets}
