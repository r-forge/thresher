\name{extractOneLGF}
\alias{RcytoGPS}
\title{extractOneLGF}
\description{
  The purpose of this function is to be used as a supplementary function for readLGF which is found within the RcytoGPS library. This function does all the parsing and formatting of the JSON file obtained from the website CytoGPS.org. extractOneLGF extracts the LGF data, along with the status of each karyotype, and generates the cytoData used for the image function from the JSON file. The status reports if the karyotype extracted was a Successâ or a Fixable grammer error and Success or if there was another error with the karyotype such as a Validation error and extracts the status of each karyotype as well as the ISCN Karyotype band. Status also uses a function called padnames to give sample row names to the statuses and karyotypes. The LGF data is extracted and converted into a data frame with new columns being generated such as the karyotype ID and Clones (rnames) along with the cytogenetic bands being produced in triple kits, each set represents the cytogenetic abnormality they are associated with. 

}
\usage{
lapply(myJSON, extractOneLGF)
}
\arguments{
  \describe{
    \item{\code{myJSON}}{ A list of elements extracted from the JSON file by using lapply}
    \item{\code{extractOneLGF}}{Essentially the FUN argument used in lapply, extractOneLGF is applied to every element of myJSON to extract and format Status and the LGF data.}
  }
}

\value{
   Returns a nested list of the binary LGF model for each element of a karyotype as well as a list of which karyotype had a 'success' or 'fail' cytogenetic event.

}
\source{

}
\references{

}
\examples{
 myJSON <- lapply(filenames, function(x) fromJSON(file = x)) # a list, one element per JSON file, uses fromJSON from the R library rjson
 temp  <- lapply(myJSON, extractOneLGF)

}
\keyword{

}