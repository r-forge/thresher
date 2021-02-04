\name{preclean}
\alias{preclean}
\title{Pre-clean Karyotypes}
\description{
  A function to clean karyotype data, by deleting known comments that do not
  adhere to the ISCN standard.
}
\usage{
preclean(x, targetColumns, dirt)
}
\arguments{
  \item{x}{A data frame containing at least one column of karyotype data.}
  \item{targetColumns}{Either a numeric vector of column indices or a
    character vector of column names.}
  \item{dirt}{A character vector containing items to delete from all karyotypes,}
}

\value{
  A data frame of te same size and with the same number of columns as
  the input data frame.
}
\details{
  The core input data worked on by the \code{RCytoGPS} are karyotypes,
  which are text strings written to conform to the ISCN standard. At
  many institutions, the cytogeneticists have developed idiosyncratic
  conventions that they use to add comments into the string. In most
  cases, these karyotypes are simply stored as text strings in a local
  database. In particular, they are not checked for synatx or grammar
  errors. By contrast, the implementation of the CytoGPS algorithm at
  the web site \url{http://cytogps.org} uses a formal approach with
  lexer and parser. As a result, many karyotypes are rejected by the
  system.

  The \code{preclean} function uses \code{\link{gsub}} to delete a list
  of known (local) comments from all karyotypes, making it more likley
  that they will be successfully processed by the lexer and parser. We
  provide an example list derived from experience at our own
  institution.

  \emph{Implementation Note:} The \code{preclean} function removes
  strings in the order that they are contained in the \code{dirt}
  vector. So, you have to be carefully not to delete parts of a long
  phrase before trying to delete the whole phrase. For example, you
  should not remove "clonal" before removing "nonclonal". 
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com}
}
\references{
  Abrams ZB, Zhang L, Abruzzo LV, Heerema NA, Li S, Dillon T, Rodriguez R,
  Coombes KR, Payne PRO. CytoGPS: a web-enabled karyotype analysis tool
  for cytogenetics. Bioinformatics. 2019 Dec 15;35(24):5365-5366.

  Abrams ZB, Tally DG,  Zhang L, Coombes CE, Payne PRO, Abruzzo LV,
  Coombes KR. Pattern recognition in lymphoid malignancies using CytoGPS
  and Mercator. In Press.

  Abrams ZB, Tally DG, Coombes KR. RCytoGPS: An R Package for Analyzing
  and Visualizing Cytogenetic Data. In preparation.
}
\examples{
cleanDir <- system.file("PreClean", package = "RCytoGPS")
bad <- read.delim(file.path(cleanDir, "badStrings.txt"), header=FALSE)
bad <- as.vector(as.matrix(bad))
input <- read.csv(file.path(cleanDir, "startKaryotypes.csv"))
myclean <- preclean(input, 4:5, bad)
}
\keyword{
character
}
