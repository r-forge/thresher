\name{cytoband-data}
\alias{cytoband-data}
\alias{cytoData}
\docType{data}
\title{Example Cytoband Data}
\description{
  This data set contains the genomic locations of cytobands along with
  loss-gain-fusion percentages for three group of samples.
}
\usage{
data(cytoData)
}
\format{A data matrix (\code{cytobandLocations}) containing 868 rows and
  14 columns. Each row contains one of the cytobands defined in the 2016
  update to ISCN nomenclature. The rownames are the standard cytoband
  names; for example, \code{1p36.33}. The first five columns are the
  same as in the \code{\link{cytobandLocations}} data set. The remaining
  nine columns form three sets of three, recording the percentage of Loss,
  Gain, and Fusion events in three sets of samples, labeled "A", "B",
  and "C".
}
\source{
  The cytoband locations were downloaded from the UCSC Genome Browser
  and synchronized with the list of cytobands in ISCN 2016. The
  percentages were computed using tools in this package from subsets of
  the Mitelman Database of Chromosome Aberrations and Gene Fusions in
  Cancer.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  Dwayne G. Tally \email{dtally110@hotmail.com}
}
\references{
  J. McGowan-Jordan, A. Simons, M. Schmid (editors).
  ISCN 2016: An International System for Human Cytogenomic Nomenclature.
  Karger Publishing, Basel, 2016.

  Mitelman, F., B. Johansson, and F. Mertens, Catalog of chromosome
  aberrations in cancer. Vol. 1. 1991: Wiley-Liss New York.
}
\keyword{datasets}

