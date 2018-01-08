\name{Thresher-class}
\Rdversion{1.1}
\docType{class}
\alias{Thresher-class}
\alias{getColors,Thresher-method}
\alias{getSplit,Thresher-method}
\alias{heat,Thresher-method}
\alias{makeFigures,Thresher-method}
\alias{plot,Thresher,missing-method}
\alias{scatter,Thresher-method}
\alias{screeplot,Thresher-method}

\title{Class \code{"Thresher"}}
\description{
The \code{Thresher} class represents the first step of an algorithm that
combines outlier detection with clustering.
}
\section{Objects from the Class}{
  Objects should be defined using the \code{\link{Thresher}} constructor.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"}; the name of
      this object. }
    \item{\code{data}:}{Object of class \code{"matrix"}; the data that
      was used for clustering. }
    \item{\code{spca}:}{Object of class \code{"SamplePCA"}; represnts
      the results of performing a principal components analysis ont he
      original \code{data}.}
    \item{\code{loadings}:}{Object of class \code{"matrix"}; the matrix
      of loading vectors from the principal components analysis.}
    \item{\code{gc}:}{Object of class \code{"hclust"}; the result of
      hierarchical clustering the data columns.}
    \item{\code{pcdim}:}{Object of class \code{"numeric"}; the number of
      significant principal components.}
    \item{\code{delta}:}{Object of class \code{"numeric"}; the lenth of
      the loading vector in the principal component space of dimension
      \code{pcdim}.}
    \item{\code{ag}:}{Object of class \code{"AuerGervini"}; represnts
      the result of running the automated Auer-Gervini algorithm to
      detemine the number of principal comonents.}
    \item{\code{agfun}:}{A function, which is used as the default method
      for computing the principal component dimension from the Auer-Gervini
      plot.}
  }
}
\section{Methods}{
  \describe{
    \item{getColors}{\code{signature(object = "Thresher")}: ... }
    \item{getSplit}{\code{signature(object = "Thresher")}: ... }
    \item{heat}{\code{signature(object = "Thresher")}: ... }
    \item{makeFigures}{\code{signature(object = "Thresher")}: ... }
    \item{plot}{\code{signature(x = "Thresher", y = "missing")}: ... }
    \item{scatter}{\code{signature(object = "Thresher")}: ... }
    \item{screeplot}{\code{signature(x = "Thresher")}: ... }
	 }
}
\references{
  Wang M, Abrams ZB, Kornblau SM, Coombes KR.
  Thresher: determining the number of clusters while removing outliers.
  BMC Bioinformatics, In press.
}
\author{
Kevin R. Coombes <krc@silicovore.com>
}
\seealso{
  \code{\link{Thresher}}, \code{\link{Reaper-class}}, \code{\link{AuerGervini-class}}
}
\examples{
showClass("Thresher")
}
\keyword{classes}
\keyword{cluster}
