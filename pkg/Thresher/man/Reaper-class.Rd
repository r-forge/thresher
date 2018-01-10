\name{Reaper-class}
\Rdversion{1.1}
\docType{class}
\alias{Reaper-class}
\alias{Reaper}
\alias{getColors,Reaper-method}
\alias{getSplit,Reaper-method}
\alias{getStyles,Reaper-method}
\alias{makeFigures,Reaper-method}

\title{Class \code{"Reaper"}}
\description{
The \code{Reaper} class implements the second step in the algorithm to
combine outlier detection with cliustering.  The first step, implemented
in the \link{Thresher-class}, performs principal components analysis an
computes the PC dimension.  Features with short loading vectors are
identified as outliers. Remaining features are clustering, based on the
directions of the loading vectors, using mixtures of von Mises-Fisher
distributions. 
}
\section{Objects from the Class}{
  Objects should be defined using the \code{Reaper} constructor. In
  the simplest case, you simply pass in a previously computed
  \code{Thresher} object.
}
\usage{
Reaper(thresher, useLoadings = FALSE, cutoff = 0.3,
       metric = NULL, linkage="ward.D2",
       maxSampleGroups = 0, ...)
}
\arguments{
  \item{thresher}{A \code{Thresher} object.}
  \item{useLoadings}{A logical value; should model-based clustering
    using von Mises-Fisher distributions be performed in the principal
    component space?}
  \item{cutoff}{A real number; what length loading vector should be used
    to separate outliers from significant contributers.}
  \item{metric}{A character string containing the name of a clustering metric
    recognized by either \code{\link{dist}} or \code{\link{distanceMatrix}}.} 
  \item{linkage}{A character string containing the name of a linkage rule
    recognized by \code{\link{hclust}}.}
  \item{maxSampleGroups}{An integer; the maximum number of sample groups
    to be indicated by color in plots of the object.}
  \item{\dots}{Additional arguments to be passed to the
    \code{\link{Thresher}} function.}
}
\section{Slots}{
  \describe{
    \item{\code{useLoadings}:}{Logical; should model-based clustering be
      performed in PC space?}
    \item{\code{keep}:}{Logical vector: which of the features (columns)
      should be retained as meaningful signal instead of being removed
      as outliers?}
    \item{\code{nGroups}:}{Object of class \code{"number or miss"}; the
      optimal number of groups/clusters found by the algorithm. If all
      of the fits fail, this is NA. }
    \item{\code{fit}:}{Object of class \code{"fit or miss"}; the best
      mixture model fit.  Can be an NA if something goes wrong
      when trying to fit mixture models.}
    \item{\code{allfits}:}{Object of class \code{"list"}; a list, each
      of whose entries should be the results of fitting a mixture model
      with a different number of components.}
    \item{\code{bic}:}{Object of class \code{"number or miss"}; the
      optimal valus of the Bayes Information Criterion; can be NA if all
      attempts to fit models fail.}
    \item{\code{metric}:}{A character string; the preferred distance
      metric for hierarchical clustering. If not specified by the user,
      then this is computed using the \code{bestMetric} function.}
    \item{\code{signalSet}:}{Object of class \code{\link{SignalSet}}}
    \item{\code{maxSampleGroups}:}{An integer; the maximum number of
      sample groups to be distinguished by color in plots of the
      object.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{Thresher}"}, directly.
}
\section{Methods}{
  \describe{
    \item{makeFigures}{\code{signature(object = "Reaper")}: This is a
      convenience function to produce a standard set of figures. In
      addition tot he plots preodcued for\code{Thresher} object, this
      function also produces heatmaps where sample clustering depends
      on either the continuous or binary signal sets.
      If the \code{DIR} argument is
      non-null, it is treated as the name of an existing directory where the
      figures are stored as PNG files. Otherwise, the figures are
      displayed interactively, one at a time, in a window on screen.}
    \item{getColors}{\code{signature(object = "Reaper")}: Returns the
      vector of colors assigned to the clustered columns in the data set.}
    \item{getSplit}{\code{signature(object = "Reaper")}:  Returns the
      vector of colors assigned to the clustered rows in the data set.}
	 }
}
\details{
  Using the dimension computed when constructing the
  \code{\link{Thresher}} object, we computed the lengths of the loading
  vectors associated to features in the data set.  Features whose length
  is less than a specified \code{cutoff} are identified as outliers and
  removed. (Based on extensive simulations, the default cutoff is
  taken to be 0.3.) We then refit the Thresher model on the remaining
  features, which should, in theory, leave the PC dimension, \code{D},
  unchanged.  We then rescale the remaining loading vectors to unit
  length, so they can be viewed as points on a hypersphere. In order to
  cluster points on a hypersphere, we use a model based on a mixture
  of von Mises-Fisher distributions. We fit mixtures for every integer
  in the range \eqn{D \le N \le 2D+1}{D <= N <= 2D+1}; this range accounts for the
  possibility that each axis has both positively and negatively
  correlated features.  The extra \eqn{+1}{+1} handles the degenerate case when
  \eqn{D=0}{D = 0}.  The best fit is determined using the Bayes Information
  Criterion (BIC). The final step is to compute a
  \code{\link{SignalSet}}; see the description of that class for more
  details.
}
\value{
  The \code{Reaper} function returns an object of the Reaper class.
}
\references{
  Wang M, Abrams ZB, Kornblau SM, Coombes KR.
  Thresher: determining the number of clusters while removing outliers.
  BMC Bioinformatics, 2018; 19(1):1-9. doi://10.1186/s12859-017-1998-9.

  Wang M, Kornblau SM, Coombes KR.
  Decomposing the Apoptosis Pathway Into Biologically Interpretable
  Principal Components.
  bioRxiv, 2017. doi://10.1101/237883.

  Banerjee A, Dhillon IS, Ghosh J, Sra S.
  Clustering on the unit hypersphere using von Mises-Fisher distributions.
  Journal of Machine Learning Research, 2005; 6:1345--1382.

  Kurt Hornik and Bettina Gr\"un.
  movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions.
  Journal of Statistical Software, 2014; 58(10):1--31.
}
\author{
 Kevin R. Coombes <krc@silicovore.com>, Min Wang.
}
\seealso{
  \code{\link{PCDimension}}, \code{\link{SignalSet}}.
}
\examples{
# Simulate  a data set with some structure
set.seed(250264)
sigma1 <- matrix(0, ncol=16, nrow=16)
sigma1[1:7, 1:7] <- 0.7
sigma1[8:14, 8:14] <- 0.3
diag(sigma1) <- 1
st <- SimThresher(sigma1, nSample=300)
# Threshing is completed; now we can reap
reap <- Reaper(st)
screeplot(reap, col='pink', lcol='red')
scatter(reap)
plot(reap)
heat(reap)
}
\keyword{classes}
\keyword{cluster}
