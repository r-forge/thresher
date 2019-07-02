\name{Mercator-class}
\docType{class}
\alias{Mercator}
\alias{addVisualization}
\title{
The \code{Mercator} Distance Visualization Object
}
\description{
The \code{Mercator} object applies one of 10 metrics of 
\link{\code{binaryDistance}} to an object of class 
\link{\code{BinaryMatrix}} and implements 5 visualizations for
clusters of large-scale, multi-dimensional data: 
hierarchical clustering, heat map, multi-dimensional scaling,
t-Stochastic Neighbor Embedding (t-SNE), and iGraph.
}
\usage{
Mercator(binaryMat, metric, method, K, ...)
addVisualization(DV, method, ...)
}
\arguments{
  \item{binaryMat}{A \link{\code{BinaryMatrix}}.}
  \item{metric}{A \link{\code{binaryDistance}} limited to the names of 
    10 selected distance metrics: \code{jaccard}, \code{sokalMichener}, \code{hamming},
    \code{russellRao}, \code{pearson}, \code{goodmanKruskal}, \code{manhattan},
    \code{canberra}, \code{binary}, or \code{euclid}.}
  \item{method}{A visualization method, limited to \code{hclust}, \code{heat}, 
    \code{mds}, \code{tsne}, and \code{graph}.}
  \item{K}{An \code{integer} specifying the number of desired clusters.}
}
\section{Slots}{
  \describe{
    \item{\code{metric}:}{Object of class \code{"character"}; the name 
      of the \link{\code{binaryDistance}} applied to create this object. }
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
  The \code{Mercator} function constructs and returns an object of 
  \code{Mercator} class, including a distance matrix calculated on a 
  given metric and a given visualiations. The \code{addVisualizations} 
  function can be used to add additional visualizations to a given 
  \code{Mercator} object.
}
\author{Kevin R. Coombes <krc@silicovore.com}
\examples{






}
\keyword{ cluster }
