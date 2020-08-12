\name{setClusters}
\alias{setClusters}
\alias{recolor}
\alias{remapColors}
\alias{recluster}
\title{
Changing Cluster and Color Assignments in Mercator
}
\description{
  Cluster assignments from unsupervised analyses typically assign
  arbitrary integers to the classes. When comparing the results of
  different algorithms or different distance metrics, it is helpful to
  match the integers in order to use colors and symbols that are as
  consistent as possible. These functions help achieve that goal.
}
\usage{
setClusters(DV, clusters)
recolor(DV, clusters)
remapColors(fix, vary)
recluster(DV, K)
}
\arguments{
  \item{DV}{An object of the \code{Mercator} class.}
  \item{clusters}{An integer vector specifiny the cluster membership of
    each element.}
  \item{fix}{An object of the \code{Mercator} class, used as the source
    of color and cluster assignments.}
  \item{vary}{An object of the \code{Mercator} class, used as the target
    of color and cluster assignments.}
  \item{K}{An integer, the number of clusters.}
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

  The functions \code{setClusters} and \code{remapColors} solve this problem
  in the context of \code{Mercator} objects. They accaomplish this task
  using the greedy algorithm implemented and described in the
  \code{\link{remap}} function in the \code{\link{Thresher}} package.
}
\value{
Both \code{setClusters} and \code{remapColors} return an object of class
\code{Mercator} in which only the cluster labels and associated color
and symbol representations (but not the distance metric used, the number
of clusters, the cluster assignments, nor the views) have been updated.

The \code{recolor} function is currently an alias for
\code{setClusters}. However, using \code{recolor} is deprecated, and the
alias will likely be removed in the next version of the package.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\seealso{ \code{\link{remap}} }
\examples{
#Form a BinaryMatrix
data("iris")
my.data <- t(as.matrix(iris[,c(1:4)]))
colnames(my.data) <- 1:ncol(my.data)
my.binmat <- BinaryMatrix(my.data)

# Form a Mercator object; Set K to the number of known species
my.vis <- Mercator(my.binmat, "euclid", "tsne", K=3)
table(getClusters(my.vis), iris$Species)
summary(my.vis)

# Recolor the Mercator object with known species
DS <- recolor(my.vis, as.numeric(iris$Species))
table(getClusters(DS), iris$Species)

# Use a different metric
my.vis2 <- Mercator(my.binmat, "manhattan", "tsne", K=3)
table(getClusters(my.vis2), iris$Species)
table(Pearson = getClusters(my.vis2), Euclid = getClusters(my.vis))

# remap colors so the two methods match as well as possible
my.vis2 <- remapColors(my.vis, my.vis2)
table(Pearson = getClusters(my.vis2), Euclid = getClusters(my.vis))

# recluster with K=4
my.vis3 <- recluster(my.vis, K = 4)

# view the results
opar <- par(mfrow=c(1,2))
plot(my.vis, view = "tsne", main="t-SNE plot, Euclid")
plot(my.vis2, view = "tsne",  main="t-SNE plot, Pearson")
par(opar)
}
\keyword{ cluster }
