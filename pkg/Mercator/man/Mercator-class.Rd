\name{Mercator-class}
\docType{class}
\alias{Mercator-class}
\alias{Mercator}
\alias{[,Mercator-method}
\alias{dim,Mercator-method}
\alias{hist,Mercator-method}
\alias{summary,Mercator-method}
\alias{plot,Mercator,missing-method}
\alias{addVisualization}
\alias{getClusters}
\title{
The \code{Mercator} Distance Visualization Object
}
\description{
  The \code{Mercator} object represents a distance matrix together with
  clustering assignments and a set of visualizations. It implements four
  visualizations for clusters of large-scale, multi-dimensional data:
  hierarchical clustering, multi-dimensional scaling, t-Stochastic
  Neighbor Embedding (t-SNE), and iGraph. The default
  \code{Mercator} constructor applies one of ten metrics of
  \code{\link{binaryDistance}} to an object of the
  \code{\link{BinaryMatrix}} class.
}
\usage{
Mercator(binaryMat, metric, method, K, ...)
addVisualization(DV, method, ...)
getClusters(DV)
}
\arguments{
  \item{binaryMat}{A \code{\link{BinaryMatrix}}.}
  \item{metric}{A \code{\link{binaryDistance}} currently limited to the
    names of 10 selected distance metrics: \code{jaccard},
    \code{sokalMichener}, \code{hamming},
    \code{russellRao}, \code{pearson}, \code{goodmanKruskal},
    \code{manhattan},
    \code{canberra}, \code{binary}, or \code{euclid}.}
  \item{method}{A visualization method, currently limited to
    \code{hclust},
    \code{mds}, \code{tsne}, and \code{graph}.}
  \item{K}{An \code{integer} specifying the number of desired clusters.}
  \item{DV}{A distance visualization produced as the output of the 
    \code{Mercator} function.}
  \item{\dots}{Additional arguments passed on to the functions that
    implement different methods for \code{addVisualization} (possibly
    passed here through the \code{Mercator} function). These include
    \describe{
      \item{\code{any}}{Any arguments to the \code{\link{cmdscale}} function for an
	\code{mds} visualization. Most commonly, this is likely to
	include \code{k}, the number of dimensions to compute.}
      \item{\code{any}}{Any arguments to the \code{\link{Rtsne}} function for a
	\code{tsne} visualization. Common examples include \code{dims}
	and \code{perplexity}.}
      \item{\code{Q}}{A quantile cutoff for the creation of the IGraph visualization.
	By default, the value is set at the 10th percentile.}
    }
  }
}
\section{Slots}{
  \describe{
    \item{\code{metric}:}{Object of class \code{"character"}; the name 
      of the \code{\link{binaryDistance}} applied to create this object. }
    \item{\code{distance}:}{Object of class \code{"dist"}; the distance 
      matrix used and represented by this object.}
    \item{\code{view}:}{Object of class \code{"list"}; contains the
      results of calculations to generate each visualize the object.}
    \item{\code{colv}:}{Object of class \code{"character"}; a vector 
      identifying the color assigned to each clustered element.}
    \item{\code{symv}:}{Object of class \code{"numeric"}; a vector 
      identifying the plotting symbol used to display each element.}
  }
}
\section{Methods}{
  \describe{
    \item{plot(x, view = NULL, ask = NULL, ...):}{
      Produce a plot of one or more visualizations within a 
      Mercator object. The default \code{view}, when omitted, is the
      first one contained in the object. You can request multiple views
      at once; if the current plot layout doesn't have enough space in
      an interactive session, the \code{ask} parameters detemines whether
      the system will ask you before advancing to the next plot. When
      plotting a graph view, you can use an optional \code{layout}
      parameter to select a specific layout by name.
    }
    \item{hist:}{\code{signature(object = "Mercator")}}{
      Produce a histogram of distances calculated in the dissimilarity 
      matrix generated in the \code{Mercator} object.
    }
    \item{summary:}{\code{signature(object = "Mercator")}}{
      Returns the chosen distance metric, dimensions of the distance 
      matrix, and available, calculated visualizations in this object.
    }
    \item{dim:}{\code{signature(object = "Mercator")}}{
      Returns the dimensions of the distance matrix of this object.
    }
    \item{[:}{\code{signature(object = "Mercator")}}{
      Subsets the distance matrix of this object.
    }
  }
}
\value{
  The \code{Mercator} function constructs and returns a distance
  visualization object of the
  \code{Mercator} class, including a distance matrix calculated on a 
  given metric and given visualizations. It is also possible (though not
  advisable) to construct a \code{Mercator} obecject directly using the
  \code{new} function.
  
  The \code{addVisualizations} function can be used to add additional 
  visualizations to an existing \code{Mercator} object.
  
  The \code{getClusters} function returns a vector of cluster assignments.
}
\author{Kevin R. Coombes <krc@silicovore.com>, Caitlin E. Coombes}
\examples{
#Form a BinaryMatrix
data("iris")
my.data <- as.matrix(iris[,c(1:4)])
my.rows <- as.data.frame(c(1:length(my.data[,1])))
my.binmat <- BinaryMatrix(my.data, , my.rows)
my.binmat <- t(my.binmat)
summary(my.binmat)

# Form a Mercator object
# Set K to the known number of species in the dataset
my.vis <- Mercator(my.binmat, "euclid", "hclust", K=3) 
summary(my.vis)
hist(my.vis)
my.vis <- addVisualization(my.vis, "mds")
plot(my.vis, view = "hclust")
plot(my.vis, view = "mds")

#Recover cluster identities
#What species comprise cluster 1?
my.clust <- getClusters(my.vis)
my.species <- iris$Species[my.clust == 1]
my.species
}
\keyword{ cluster }
