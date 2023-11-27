\name{readLGF}
\alias{readLGF}
\title{Extracting LGF karyotype data from JSON files }
\description{
  A function to read binary karyotype data, stored in LGF format in JSON
  files produced by the CytoGPS web site, into R for further analysis.
}
\usage{
readLGF(files = NULL, folder = NULL, verbose = TRUE)
}
\arguments{
  \item{files}{The name of the JSON file (or a character vector of
    such file names) from which you want to extract and format
    data. If \code{NULL}, then it will extract all JSON files within
    the folder path provided.}
  \item{folder}{The specified directory/folder from which the user
    wants to extract JSON files. If \code{NULL}, then the function
    will look in the current working directory.}
  \item{verbose}{A logical value; should the function keep you informed about what it is doing?}
}

\value{
  A list containing five elements:

    \item{\code{source}:}{A character vector containing the names of the
      JSON files from which data was read.}
    \item{\code{raw}:}{A list of lists, one per JSON source file. Each
      internal list contains two elements, \code{Status} and
      \code{LGF}. \code{Status} is a data frame with one row per
      karyotype in the input file, describing the results of CytoGPS
      parsing and mapping. Results can be "\code{Success}",
      "\code{Nonfixable grammar error}",
      "\code{Validation error}",
      "\code{Fixable grammar error and success}",  or
      "\code{Fixable grammar error but containing validation error}".
      \code{LGF} is a data frame where the columns are LGF-cytobands and the
      rows are clones from the successfully processed input karyotypes;
      each karyotype can have multiple clones. Entries are zero or one
      indicating the absence or presence of an abnormality.}
    \item{\code{frequency}:}{A data frame, where the rows are cytobands
      and the columns are Loss, Gain, and Fusion, with three columns per
      input file. Entries are the fraction of karyotype clones with
      that abnormality.}
    \item{\code{size}:}{An integer vector containing the total number of
      clones per input file. These values can be used to turn fractions
      back into counts.}
    \item{\code{CL}:}{A data frame with one row per cytoband detailing
      the chromosomal location and (grayscale) color of the band
      produced by Giemsa staining.}

}
\details{

  CytoGPS is an algorithm that converts conventional karyotypes from the
  standard text-based notation (the International Standard for Human
  Cytogenetic/Cytogenomic Nomenclature; ISCN) into binary vectors with
  three bits (loss, gain, or fusion) per cytoband, which we call the LGF
  model. The web site \url{http://cytogps.org} provides an
  implementation that allows users to upload text files containing one
  karyotype per line. It produces its output as a file in JavaScript
  Object Notation (JSON).

  The \code{readLGF} function reads and parses these JSON files and
  converts them into an R data structure. The \code{raw} component of
  this structure contains binary matrices that can serve as input to the
  Mercator package (see \code{\link[Mercator]{Mercator-class}}) for
  unsupervised analyses. The \code{frequency} component summarizes the
  fraction of input karyotype-clones with each abnormality, and can be
  visualized with other tools in the \code{RCytoGPS} package.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  Dwayne G. Tally \email{dtally110@hotmail.com}
}
\seealso{
  \code{\link[Mercator]{Mercator-class}}
}
\references{
  Abrams ZB, Zhang L, Abruzzo LV, Heerema NA, Li S, Dillon T, Rodriguez R,
  Coombes KR, Payne PRO. CytoGPS: a web-enabled karyotype analysis tool
  for cytogenetics. Bioinformatics. 2019 Dec 15;35(24):5365-5366.

  Abrams ZB, Tally DG,  Zhang L, Coombes CE, Payne PRO, Abruzzo LV,
  Coombes KR. Pattern recognition in lymphoid malignancies using CytoGPS
  and Mercator. Under review.

  Abrams ZB, Tally DG, Coombes KR. RCytoGPS: An R Package for Analyzing
  and Visualizing Cytogenetic Data. In preparation.
}
\examples{
jsonDir <- system.file("Examples/JSONfiles", package = "RCytoGPS")
x <- readLGF(folder = jsonDir)

jsonFile <- dir(jsonDir, pattern = "*.json")[1]
y <- readLGF(jsonFile, jsonDir)

}
\keyword{
IO
}
