\name{CytobandData-class}
\alias{CytobandData-class}
\alias{CytobandData}
\alias{summary,CytobandData-method}
\title{The \code{CytobandData} Class}
\description{
  A class to represent genome-wide data measured or summarized at
  cytoband-level resoluton.
}
\usage{
CytobandData(data, info, genome = NULL)
\S4method{summary}{CytobandData}(object, \dots)
}
\arguments{
  \item{data}{ A data frame that must contain at least one
    numeric column and may contain the five cytoband location
    columns. If the latter are missing, they must be supplied by the
    \code{genome} argument.}
  \item{info}{ A data frame with \code{Label} and
    \code{Description} columns. If missing, it will be created from the
    column names in the \code{data} argument.}
  \item{genome}{ A data frame containing cytoband
    locations. Often, simply uses the \code{\link{cytobandLocations}}
    object describing build 38 of the human genome.}
  \item{object}{An object of the \code{CytobandData} class.}
  \item{\dots}{Additional arguments to the summary method; ignored.}
}
\section{Slots}{
  \describe{
    \item{\code{DATA}:}{A data frame containg five columns  (Chromosome,
      loc.start, loc.end, Band, and Stain) defining chromosomal
      locations of cytobands, along with one or more columns of
      numerical data at cytoband resolution.}
    \item{\code{INFO}:}{A data frame with two columns (Label and
      Description) that describes the columns of the \code{DATA} slot.}
  }
}
\section{Methods}{
  \describe{
    \item{summary:}{\code{signature(object = "RCytoGPS")}}{
      Returns a summary of the \code{DATA} slot.
    }
  }
}
\details{
  Conventional karyotyping describes chromosomal abnormalities (in a
  standardized text-based nomenclature) that are visible through a
  microsope. Karyotyping relies on a technique known as Giemsa staining,
  which creates a banding pattern along the chromosome of different
  shades of gray. The Loss-Gain-Fusion (LGF) model implemented at the
  \url{http://cytogps.org} web site converts text-based karyotypes into
  binary vectors, stored in JSON files. Our \code{\link{readLGF}}
  function parses the JSON files to create R data structures, including
  cytoband-resolution summaries of the frequencies of
  abnormalities. These summaries can be used to create
  \code{CytobandData} objects, which can be visualized using the
  \code{\link{barplot}} and \code{\link{image}} methods.
}
\value{
 The \code{CytobandData} conxtructor returns an object of the
 \code{CytobandData} class. The \code{summary} method returns a table,
 which is the usual return value when the summary method is applied to a
 data frame (in this case, the \code{DATA} slot).

}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  Dwayne G. Tally \email{dtally110@hotmail.com}
}
\references{
  Abrams ZB, Tally DG, Coombes KR. RCytoGPS: An R Package for Analyzing
  and Visualizing Cytogenetic Data. In preparation.

  Abrams ZB, Tally DG,  Zhang L, Coombes CE, Payne PRO, Abruzzo LV,
  Coombes KR. Pattern recognition in lymphoid malignancies using CytoGPS
  and Mercator. Under review.
}
\examples{
jsonDir <- system.file("Examples/JSONfiles", package = "RCytoGPS")
temp <- readLGF(folder = jsonDir)
cytoData <- data.frame(temp[["CL"]],temp[["frequency"]])
bandData <- CytobandData(cytoData)
class(bandData)
summary(bandData)
}
