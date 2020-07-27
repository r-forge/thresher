\name{readLGF}
\alias{RcytoGPS}
\alias{Extracting JSON for LGF model}
\alias{RcytoGPS}
\docType{data}
\title{readLGF}
\description{
 A function that converts conventional karyotypes from the standard text-based notation (the International Standard for Human Cytogenetic Nomenclature; ISCN) into a binary vector with three bits (loss, gain, or fusion) per cytoband, which we call the LGF model. There are several supplementary functions that readLGF uses called extractOneLGF, rnames, which formats the LGF model and generates the IDS, and Idioformat; which obtains the frequency of cytogenetic events. ReadLGF also uses a class called LGF to format the output into a list that is similiar to a S4 object.
}
\usage{
readLGF(files = x, folder = y, verbose = TRUE)
}
\arguments{
  \describe{
    \item{\code{files}}{The name of the JSON file you want to extract and format data from. If not specified then it will extract all JSON files within the folder path provided.}
    \item{\code{folder}}{The specified directory/folder the user wants the function to extract JSON files from. If an absolute path was not specified, then the function will look in the current working directory.}
  }
}

\value{
     A nested list that has source, raw, freqeuncy, size, and CL organized in a s4 class.

}
\source{

}
\references{

}
\examples{
 y <- readLGF("results.JSON", wd)
 x <- readLGF(wd)

}
\keyword{

}