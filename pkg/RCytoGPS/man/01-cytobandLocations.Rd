\name{cytobandLocations}
\alias{cytobandLocations}
\alias{idiocolors}
\docType{data}
\title{Cytoband Locations}
\description{
  This data set contains the genomic locations of cytobands based on
  both the latest build of the human genome (GRch38) and the latest
  update to the International Standard for human Cytogenomic
  Nomenclature (ISCN). Note that the CRch38 locations are unchanged from
  build GRch37 (hg19). 
}
\usage{
data(cytobandLocations)
}
\format{A data matrix (\code{cytobandLocations}) containing 868 rows and
  5 columns. Each row contains one of the cytobands defined in the 2016
  update to ISCN nomenclature. The rownames are the standard cytoband
  names; for example, \code{1p36.33}. The columns are
  \describe{
    \item{Chromosome}{The name of the human chromosome, stored as
      \code{chr#}.}
    \item{\code{loc.start}}{The starting base position of the band.}
    \item{\code{loc.end}}{The ending base position of the band.}
    \item{\code{Band}}{The band name without the chromosome; for
      example, \code{p36.33}.}
    \item{\code{Stain}}{A factor containing the name of the Giemsa-stain
      color of the band in a karyotype image.}
  }
  Also, a vector (\code{idiocolors}) of length eight that translates the
  Giemsa stain names into appropriate colors.
}
\source{
  The starting point for these data were the cytoband locations
  downloaded from the UCSC Genome Browser. We confirmed that the data
  were unchanged between human genome builds GRCh36 (hg18), GRCh37
  (hg19), and GRCh38. However, the list of cytobands from UCSC was not
  consistent with the list of cytobands in ISCN 2016. We manually edited
  the source file to be compliant. It also matches the list of cytobands
  produced at \url{http://cytogps.org}.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  Dwayne G. Tally \email{dtally110@hotmail.com}
}
\references{
  J. McGowan-Jordan, A. Simons, M. Schmid (editors).
  ISCN 2016: An International System for Human Cytogenomic Nomenclature.
  Karger Publishing, Basel, 2016.
}
\keyword{datasets}

