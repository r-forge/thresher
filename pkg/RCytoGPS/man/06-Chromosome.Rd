\name{Chromosome-class}
\alias{Chromosome-class}
\alias{Chromosome}
\alias{image,Chromosome-method}
\title{The \code{Chromosome} Class}
\description{
  A class to represent a single chromsome in order to plot an image of
  the Giemsa-stained cytobands.
}
\usage{
Chromosome(I, res = 2500, maxbase = 250000000)
\S4method{image}{Chromosome}(x, horiz = FALSE, mai = NULL, showBandNames = FALSE, ...)
}
\arguments{
  \item{I}{A human chromosome identifier; one of the values in \code{c(1:22,"X", "Y")}.}
  \item{res}{An integer resolution; the number of pixels used to
    represent an entire chromosome.}
  \item{maxbase}{An integer representing an upper bound on the length of
    the longest chromosome, mneasured in base-pairs.}
  \item{x}{An object of the \code{Chromosome} class.}
  \item{horiz}{A logical value: should the image of the chromsome be
    oriented horizontally.}
  \item{mai}{margin inches, as in the usual graphical argument to \code{par}.}
  \item{showBandNames}{logical; shold the names of the bands be written on the plot?}
  \item{\dots}{Additional arguments to the image method; ignored.}
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{A character value, typically of the form \code{"chr1"}.}
    \item{\code{label}:}{A chartacter value, typiclly of the form
      \code{"Chr 1"}.}
    \item{\code{grid}:}{An integer vector (of length \code{res})
      repesenting the base position along the chromosome.}
    \item{\code{range}:}{An integer vector of length 2 marking teh
      starting and ending position of the banded part of the chromosome,
      in bases.}
    \item{\code{stain}:}{An integer vector (of length \code{res}), where
      the values indicate the color of the Giemsa stain for that part of
      the chromosome.}
  }
}
\section{Methods}{
  \describe{
    \item{image:}{\code{signature(object = "RCytoGPS")}}{
      Creates an image of the chromosme, with bands colored according
      to Giemsa staining.
    }
  }
}
\details{
  Conventional karyotyping describes chromosomal abnormalities (in a
  standardized text-based nomenclature) that are visible through a
  microsope. Karyotyping relies on a technique known as Giemsa staining,
  which creates a banding pattern along the chromosome of different
  shades of gray. This class is used to plot images of individual
  chromosomes, corectly reflecting the size and color of the bands..
}
\value{
 The \code{Chromosome} constructor returns an object of the
 \code{Chromosome} class. The \code{image} method invisibly returns its
 first argument.
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
x <- Chromosome(6)
image(x)
image(x, showBandNames = TRUE)
image(x, horiz = TRUE)
}
