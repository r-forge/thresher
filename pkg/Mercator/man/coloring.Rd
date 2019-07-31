\name{recolor}
\alias{recolor}
\alias{remapColors}
\title{
Changing Cluster and Color Assignments in Mercator
}
\description{
  Cluster assignments from unsupervised analyses typically assign
  arbitrary integers to the classes. When comparing the results of
  different algorithms or different distance metrics, it is helpful to
  match the integers in order to use colors and symbois that are as
  consistent as possible. These functions help achieve that goal.
}
\usage{
recolor(DV, clusters)
remapColors(fix, vary)
}
\arguments{
  \item{DV}{An object of eh \code{Mercator} class.}
  \item{clusters}{An integer vector specifiny the cluster membership of
    each element.}
  \item{fix}{An object of the \code{Mercator} class, used as the source
    of color and cluster assignments.}
  \item{vary}{An object of the \code{Mercator} class, used as the target
    of color and cluster assignments.}
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

  The functions \code{recolor} and \code{remapColors} solve this problem
  in the context of \code{Mercator} objects. They accaomplish this task
  using the greedy algorithm implemented and described in the
  \code{\link{remap}} function in the \code{\link{Thresher}} package.
}
\value{
Both \code{recolor} and \code{remapColors} return an object of class
\code{Mercator} in which only the cluster labels and associated color
and symbol representations (but not the distance metric used, the number
of clusters, the cluster assignmewnts, nor the views) have been updated.
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

# view the results
opar <- par(mfrow=c(1,2))
plot(my.vis@view[[1]]$Y, col=my.vis@colv, pch=my.vis@symv,
     main="t-SNE plot, Euclid",
     xlab="T1", ylab="T2")
plot(my.vis2@view[[1]]$Y, col=my.vis2@colv, pch=my.vis2@symv,
     main="t-SNE plot, Pearson",
     xlab="T1", ylab="T2")
par(opar)
}
\keyword{ cluster }
