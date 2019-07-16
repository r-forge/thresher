\name{Mercator-class}
\docType{class}
\alias{Mercator}
\alias{addVisualization}
\title{
The \code{Mercator} Distance Visualization Object
}
\description{
The \code{Mercator} object applies one of 10 metrics of 
\code{\link{binaryDistance}} to an object of class 
\code{\link{BinaryMatrix}} and implements 5 visualizations for
clusters of large-scale, multi-dimensional data: 
hierarchical clustering, heat map, multi-dimensional scaling,
t-Stochastic Neighbor Embedding (t-SNE), and iGraph.
}
\usage{
Mercator(binaryMat, metric, method, K, ...)
addVisualization(DV, method, ...)
getClusters(DV)
}
\arguments{
  \item{binaryMat}{A \code{\link{BinaryMatrix}}.}
  \item{metric}{A \code{\link{binaryDistance}} limited to the names of 
    10 selected distance metrics: \code{jaccard}, \code{sokalMichener}, \code{hamming},
    \code{russellRao}, \code{pearson}, \code{goodmanKruskal}, \code{manhattan},
    \code{canberra}, \code{binary}, or \code{euclid}.}
  \item{method}{A visualization method, limited to \code{hclust}, \code{heat}, 
    \code{mds}, \code{tsne}, and \code{graph}.}
  \item{K}{An \code{integer} specifying the number of desired clusters.}
  \item{DV}{A distance visualization produced as the output of the 
    \code{Mercator} function.}
  \item{Q}{A quantile cutoff for the creation of the IGraph visualization. 
    By default, the value is set at the 10th percentile.}
  }
}
\section{Slots}{
  \describe{
    \item{\code{metric}:}{Object of class \code{"character"}; the name 
      of the \code{\link{binaryDistance}} applied to create this object. }
    \item{\code{distance}:}{Object of class \code{"dist"}; the distance 
      matrix calculated by this object.}
    \item{\code{view}:}{Object of class \code{"list"}; contains the 
      calculations to generate each visualization added to the object.}
    \item{\code{colv}:}{Object of class \code{"character"}; a vector 
      identifying the color assigned to each clustered feature.}
    \item{\code{symv}:}{Object of class \code{"numeric"}; a vector 
      identifying the cluster number assigned to each clustered feature.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(object = "Mercator"@view[[n]], col = "Mercator"@colv, pch = "Mercator"@symv)}}:{
      Produce a plot of the \code{nth} visualization created within a 
      Mercator object. Producing a t-SNE plot requires the added 
      information \code{"Mercator"@view[[n]]$Y}
    }
    \item{hist}{\code{signature(object = "Mercator")}}:{
      Produce a histogram of distances calculated in the dissimilarity 
      matrix generated in the \code{Mercator} object.
    }
    \item{summary}{\code{signature(object = "Mercator")}}:{
      Returns the chosen distance metric, dimensions of the distance 
      matrix, and available, calculated visualizations in this object.
    }
    \item{dim}{\code{signature(object = "Mercator")}}:{
      Returns the dimensions of the distance matrix of this object.
    }
    \item{[}{\code{signature(object = "Mercator")}}:{
      Subsets the distance matrix of this object.
    }
  }
}
\value{
  The \code{Mercator} function constructs and returns a distance
  visualization object of 
  \code{Mercator} class, including a distance matrix calculated on a 
  given metric and given visualiations.
  
  The \code{addVisualizations} function can be used to add additional 
  visualizations to a given \code{Mercator} object.
  
  The \code{getClusters} function returns the cluster number and color 
  assigned to each observation in the \code{Mercator} distance 
  visualization.
}
\author{Kevin R. Coombes <krc@silicovore.com}
\examples{
#Form a BinaryMatrix
data("iris")
my.data <- as.matrix(iris[,c(1:4)])
my.rows <- as.data.frame(c(1:length(my.data[,1])))
my.binmat <- BinaryMatrix(my.data, , my.rows)
my.binmat <- t(my.binmat)
summary(my.binmat)

#Form a Mercator object
my.vis <- Mercator(my.binmat, "euclid", "hclust", K=3) #set K to the known number of species in the dataset
summary(my.vis)
hist(my.vis)
my.vis <- addVisualization(my.vis, "mds")
plot(my.vis@view[[1]], col=my.vis@colv, pch=my.vis@symv)
plot(my.vis@view[[2]], col=my.vis@colv, pch=my.vis@symv)

#Recover cluster identities and colors
my.clust <- getClusters(my.vis)
my.clust






}
\keyword{ cluster }
