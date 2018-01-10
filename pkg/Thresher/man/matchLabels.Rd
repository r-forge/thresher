\name{matchLabels}
\alias{matchLabels}
\alias{labelMatcher}
\alias{labelAccuracy}
\alias{countAgreement}
\alias{bestMetric}
\title{
Match Arbitrary Class Assignments Across Methods
}
\description{
These functions provide a set of tools to find the best match between
the labels used by two different algorithms to cluster the same set of
samples.
}
\usage{
labelMatcher(tab, verbose = FALSE)
matchLabels(tab)
countAgreement(tab)
labelAccuracy(data, labels, linkage="ward.D2")
bestMetric(data, labels)
}
\arguments{
  \item{tab}{A contingency table, represented as a square \code{matrix}
    or \code{table} as an R object.  Both dimensions represent an
    assignment of class labels, with each row and column representing
    one of the labels.  Entries should be non-negative integer counts of
    the number of objects having the labels represented by the row and
    column.}
  \item{verbose}{A logical value; should the routine print something out
    periodically so you know it's still working?}
  \item{data}{A matrix whose columns represent objects to be clustered
    and whose rows represent the anonymous features used to perform the
    clustering.}
  \item{labels}{A factor (or character vector) of class labels for the
    objects in the \code{data} matrix.}
  \item{linkage}{A linkage rule accepted by the \code{\link{hclust}}
    function.}
}
\details{
  In the most general sense, clustering can be viewed as a function from
  the space of "objects" of interest into a space of "class labels".  In
  less mathematical terms, this simply means that each object gets
  assigned an (arbitrary) class label. This is all well-and-good until
  you try to compare the results of running two different clustering
  algorithms that use different labels (or even worse, use the same
  labels -- typically the integers \eqn{1, 2, \dots, K} -- with
  different meanings).  When that happens, you need a way to decide
  which labels from the different sets are closest to meaning the
  "same thing".

  That's where this set of functions comes in. The core algorithm is
  implemented in the function \code{labelMatcher}, which works on a
  contingency table whose entries \eqn{N_{ij}} are the number of samples
  with row-label = \eqn{i} and column-label = \eqn{j}.  To find  the
  best match, one computes (heuristically) the values \eqn{F_{ij}} that
  describe the fraction of all entries in row \eqn{i} and column \eqn{j}
  represented by \eqn{N_{ij}}.  Perfectly matched labels would consist
  of a row \eqn{i} and a column \eqn{j} where \eqn{N_{ij}} is the only
  nonzero entry in its row and column, so \eqn{F_{ij} = 1}. The largest
  value for \eqn{F_{ij}} (with ties broken simply by which entry is
  closer to the upper-left corner of the matrix) defines the best
  match. The matched row and column are then removed from the matrix and
  the process repeats recursively.

  We apply this method to determine which distance metric, when used in
  hierarchical clustering, best matches a "gold standard" set of class
  labels.  (These may not really be gold, of course; they can also be a
  set of labels determined by k-means or another clustering algorithm.)
  The idea is to cluster the samples using a variety of different
  metrics, and select the one whose label assignments best macth the
  standard.
}
\value{
  The \code{labelMatcher} function returns a list of two vectors of the
  same length.  These contain the matched label-indices, in the order
  they were matched by the algorithm.

  The \code{matchLabels} function is a user-friendly front-end to the
  \code{labelmatcher} function.  It returns a matrix, with the rows and
  columns reordered so the labels match.

  The \code{countAgreement} function returns an integer, the number of
  samples with the "same" labels, computed by summing the diagonal of
  the reordered matrix produced by \code{matchLabels}.

  The \code{labelAccuracy} function returns a vector indexed by the set
  of nine distance metrics hard-coded in the function.  Each entry is
  the fraction of samples whose hierarchical clusters match the
  prespecified \code{labels}.

  The \code{bestMetric} function is a user-friendly front-end to the
  \code{labelAccuracy} function.  It returns the name of the distance
  metric whose hierarchical clusters best match the prespecified
  \code{labels}. 
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\note{The \code{labelAccuracy} function should probably allow the user
  to supply a list of distance metrics instead of relying on the
  hard-coded list internally.}
\seealso{
  Hierarchical clustering is implemented in the \code{\link{hclust}}
  function.  We use the extended set of distance metrics provided by the
  \code{\link{distanceMatrix}} function from the ClassDiscovery package.
  This set includes all of the  metrics from the \code{\link{dist}}
  funciton. 
}
\examples{
factor1 <- sample(c("A", "B", "C"), 30, replace=TRUE)
factor2 <- rep(c("X", "Y", "Z"), each=10)
tab <- table(factor1, factor2)
matchLabels(tab)
labelMatcher(tab)
}
\keyword{ cluster }
\keyword{ classif }
\keyword{ multivariate }
