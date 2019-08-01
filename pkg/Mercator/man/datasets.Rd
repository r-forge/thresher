\name{mercator-data}
\alias{mercator-data}
\alias{CML500}
\alias{CML1000}
\alias{lgfFeatures}
\docType{data}
\title{CML Cytogenetic Data}
\description{
  These data sets contain binary versions of subsets of cytogenetic
  karyotype data from patients with chronic myelogenous leukemia (CML). 
}
\usage{
data("lgfFeatures")
data("CML500")
data("CML1000")
}
\format{
  \describe{
    \item{\code{lgfFeatures}}{A data matrix with 2748 rows and 6 columns
    listing the cytogentic bands produced as output of the CytoGPS
    algorithm that converts text-based karyotypes into a binary
    Loss-Gain-Fusion (LGF) model. The columns include the \code{Label}
    (the Type and Band, joined by an underscore),
    \code{Type} (Loss, Gain, or Fusion), \code{Band} (standard name of
    the cytogenetic band), \code{Chr} (chromosome), \code{Arm} (the
    chromsome arm, of the form #p or #q), and \code{Index} (an integer
    that can be used for sorting or indexing).
  }
    \item{\code{CML500}}{A \code{\link{BinaryMatrix}} object with 770
      rows (subset of LGF features) and 511 columns (patients). The
      patients were selected using the \code{\link{downsample}} function
      from the full set of more than 3000 CML karyotypes. The rows were
      selected by removing redundant and non-informative features when
      considering the full data set.}
    \item{\code{CML1000}}{A \code{\link{BinaryMatrix}} object with 770
      rows (subset of LGF features) and 1057 columns (patients). The
      patients were selected using the \code{\link{downsample}} function
      from the full set of more than 3000 CML karyotypes. The rows were
      selected by removing redundant and non-informative features when
      considering the full data set.}
  }
} 
\source{
  The cytogenetic data were obtained from the public Mitelman Database
  of Chromosomal Aberrations and Gene Fusions in Cancer on 4 April
  2019. The database is currently located at
  https://cgap.nci.nih.gov/Chromosomes/Mitelman as part of hte
  Cancer Genome Anatomy Project (CGAP). The CGAP web site is expected to
  close on 1 October 2019 at which point the Mitelman database will
  move to an as-yet-undisclosed location. The data were then converted
  from text-based karyotrypes into binary vectors using CytoGPS
  \url{http://cytogps.org/}.
}
\references{
  Abrams ZB, Zhang L, Abruzzo LV, Heerema NA, Li S, Dillon T, Rodriguez
  R, Coombes KR, Payne PRO. CytoGPS: A Web-Enabled Karyotype Analysis
  Tool for Cytogenetics. Bioinformatics. 2019 Jul 2. pii: btz520.  doi:
  10.1093/bioinformatics/btz520. [Epub ahead of print]
}
\author{
  Kevin R. Coombes <krc@silicovore.com>, Caitlin E. Coombes
}
\keyword{datasets}

