\name{downsample}
\alias{downsample}
\title{
Downsampling a Distance Visualization to Facilitate iGraph
}
\description{
The \code{downsample} function implements a structured reduction 
of data points with a parent \code{Mercator} distance visualization to 
improve visualization and computational time, especially for the 
implementation of the iGraph visualization.
}
\usage{
downsample(target, distanceMat, cutoff)
}
\arguments{
  \item{target}{An \code{integer} number of points to which the user wishes to 
    reduce the parent \code{Mercator} object.
    }
  \item{distanceMat}{An object of class \code{matrix} containing the 
    distance matrix component of the parent \code{Mercator} object.
    }
  \item{cutoff}{An inclusion cutoff for selected points based on the 
    local density within the parent data.
    }
}
\details{
Mercator can be used to visualize complex networks using 
iGraph. To improve clarity of the visualization and computational 
time, we implement the \code{downsample} function to reduce the 
number of data points to be linked and visualized. The conceptual 
grounding for \code{downsample} draws on Peng Qiu's implementation of the
SPADE clustering algorithm for mass cytometry data. The
\code{downsample} function 
under-samples the densest regions of the data space to make it more 
likely that rarer clusters will still be adequately sampled.
}
\value{
\code{downsample} returns an object of class \code{Mercator} containing 
a structured subset of items from the parent \code{Mercator} object.
}
\references{
Qiu, P., et. al. (2011). Extracting a cellular hierarchy from high-dimensional cytometry data with SPADE. Nature biotechnology, 29(10), 886.
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

# Form and plot Mercator object
# Set K to the known number of species in the dataset
my.vis <- Mercator(my.binmat, "euclid", "tsne", K=3)
summary(my.vis)
plot(my.vis@view[[1]]$Y, col=my.vis@colv, pch=my.vis@symv,
     main="t-SNE plot of all data points",
     xlab="T1", ylab="T2")

#Downsample the Mercator object
M <- as.matrix(my.vis@distance)
set.seed(21340)
DS <- downsample(50, M, 0.1)
red.vis <- my.vis[DS]

#Visualize the down sampled t-SNE plot
plot(red.vis@view[[1]]$Y, col=red.vis@colv, pch=red.vis@symv,
     main="Down sampled t-SNE Plot",
     xlab="T1", ylab="T2")
}
\keyword{ cluster }
