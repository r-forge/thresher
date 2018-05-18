\name{Thresher-class}
\Rdversion{1.1}
\docType{class}
\alias{Thresher-class}
\alias{Thresher}
\alias{heat,Thresher-method}
\alias{makeFigures,Thresher-method}
\alias{plot,Thresher,missing-method}
\alias{scatter,Thresher-method}
\alias{screeplot,Thresher-method}
\alias{getColors,Thresher-method}
\alias{getSplit,Thresher-method}
\alias{getStyles,Thresher-method}

\title{Class \code{"Thresher"}}
\description{
The \code{Thresher} class represents the first step of an algorithm that
combines outlier detection with clustering. The object combines the
results of hierarchical clustering and principal components analysis
(with a computation of its dimension) on the same data set.
}
\section{Objects from the Class}{
  Objects should be defined using the \code{Thresher} constructor. In
  the simplest case, you simply pass in the data matrix that you want to
  cluster using the Thresher algorithm.
}
\usage{
Thresher(data, nm = deparse(substitute(data)),
         metric = "pearson", linkage="ward.D2",
         method = c("auer.gervini", "broken.stick"),
         scale = TRUE, agfun = agDimTwiceMean)
}
\arguments{
  \item{data}{A data matrix.}
  \item{nm}{A character string; the name of this object.}
  \item{metric}{A character string containing the name of a clustering metric
    recognized by either \code{\link{dist}} or
    \code{\link{distanceMatrix}}.} 
  \item{linkage}{A character string containing the name of a linkage rule
    recognized by \code{\link{hclust}}.}
  \item{method}{A character string describing the algorthim used from the
    \code{PCDimension} package to compute the number of
    significant components.}
  \item{scale}{A logical value; should the data be scaled before use?}
  \item{agfun}{A function that will be accepted by the
    \code{\link{AuerGervini}} function in the
    \code{PCDimension} package.}
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"}; the name of
      this object. }
    \item{\code{data}:}{Object of class \code{"matrix"}; the data that
      was used for clustering. }
    \item{\code{spca}:}{Object of class \code{"SamplePCA"}; represents
      the results of performing a principal components analysis on the
      original \code{data}.}
    \item{\code{loadings}:}{Object of class \code{"matrix"}; the matrix
      of loading vectors from the principal components analysis.}
    \item{\code{gc}:}{Object of class \code{"hclust"}; the result of performing
      hierarchical clustering on the data columns.}
    \item{\code{pcdim}:}{Object of class \code{"numeric"}; the number of
      significant principal components.}
    \item{\code{delta}:}{Object of class \code{"numeric"}; the lengths of
      the loading vectors in the principal component space of dimension
      equal to \code{pcdim}.}
    \item{\code{ag}:}{Object of class \code{"AuerGervini"}; represents
      the result of running the automated Auer-Gervini algorithm to
      detemine the number of principal components.}
    \item{\code{agfun}:}{A function, which is used as the default method
      for computing the principal component dimension from the Auer-Gervini
      plot.}
  }
}
\section{Methods}{
  \describe{
    \item{screeplot}{\code{signature(x = "Thresher")}: Produce a scree
      plot of the PCA part of the Thresher object.}
    \item{scatter}{\code{signature(object = "Thresher")}: Produce a
      scatter plot of the first two principal components.}
    \item{plot}{\code{signature(x = "Thresher", y = "missing")}: In two
      dimensions, plot the loading vectors of the PCA part of the
      object.}
    \item{heat}{\code{signature(object = "Thresher")}: Produce a heatmap
      of the data set.}
    \item{makeFigures}{\code{signature(object = "Thresher")}: This is a
      convenience function to produce a standard set of figures for a
      \code{Thresher} object. These are (1) a scree plot, (2) a plot of
      teh Auer-Gervini slot, (3) a scatter plot of the firtst trwo
      principal components, (4) one or more plots of the loading
      vectors, depending on the PCV dimension, and (5) a heat map.
      If the \code{DIR} argument is
      non-null, it is treated as the name of an existing directory where the
      figures are stored as PNG files. Otherwise, the figures are
      displayed interactively, one at a time, in a window on screen.}
    \item{getColors}{\code{signature(object = "Thresher")}: Returns the
      vector of colors assigned to the clustered columns in the data set.}
    \item{getSplit}{\code{signature(object = "Thresher")}:  Returns the
      vector of colors assigned to the clustered rows in the data set.}
    \item{getStyles}{\code{signature(object = "Thresher")}: I refuse to
      document this, since I am not convinced that it should actually exist.}
  }
}
\details{
  \code{Thresher} operates on a data matrix that is assumed to be
  organized with rows equal to samples and columns equal to features
  (like genes or proteins).   The algorithm begins by centering and (by
  default, though this can be overridden with the \code{scale}
  parameter) standardizes the data columns. It then performs a principal
  components analysis, and uses the Auer-Gervini method, as automated in
  the \code{\link{PCDimension}} package, to determine the number, \code{D}, of
  statistically significant principal components. For each
  column-feature, it computes and remembers the length of its loading
  vector in D-dimensional space. (In case the Auer-Gervini method finds
  that \eqn{D=0}{D = 0}, the length is instead computed using \eqn{D=1}{D
    = 1}.) These
  loading-lengths will be used later to identify and remove features
  that act as outliers and do not contribute to clustering the
  samples. Finally, \code{Thresher} computes and saves the results of
  hierarchically clustering the features in the data set, using the
  specified distance \code{metric} and \code{linkage} rule.
}
\value{
  The \code{Thresher} function constructs and returns an object of the
  \code{Thresher} class.
}
\references{
  Wang M, Abrams ZB, Kornblau SM, Coombes KR.
  Thresher: determining the number of clusters while removing outliers.
  BMC Bioinformatics, 2018; 19(1):1-9. doi://10.1186/s12859-017-1998-9.

  Wang M, Kornblau SM, Coombes KR.
  Decomposing the Apoptosis Pathway Into Biologically Interpretable
  Principal Components.
  bioRxiv, 2017. doi://10.1101/237883.
}
\author{
Kevin R. Coombes <krc@silicovore.com>, Min Wang.
}
\seealso{
  \code{\link{Thresher}}, \code{\link{Reaper-class}}, \code{\link{AuerGervini-class}}
}
\examples{
set.seed(3928270)
ranData <- matrix(rnorm(100*12), ncol=12)
colnames(ranData) <- paste("G", 1:12, sep='')
thresh <- Thresher(ranData) # fit the model
screeplot(thresh)           # check the scree plot; suggests dim = 4
plot(thresh@ag, list(thresh@agfun)) # Auer-Gervini object; dim = 0
scatter(thresh)             # PCA scatter plot  (rows = samples)
plot(thresh)                # PCA loadings plot (cols = features)
heat(thresh)                # ubiquitous 2-way heatmap
}
\keyword{classes}
\keyword{cluster}
