\name{BinaryMatrix-class}
\doctype{class}
\alias{BinaryMatrix}
\alias{BinaryMatrix-class}
\alias{validBinaryMatrix}
\alias{removeDuplicateFeatures}

\title{Class \code{"BinaryMatrix"}}
\description{
The \code{BinaryMatrix} object class underlies the \code{threshLGF} and \code{Mercator} 
methods and visualizations of the \code{Mercator} package. The \code{BinaryMatrix} function
returns a new object of \code{BinaryMatrix} class.
}
\section{Objects from the Class}{
  Objects should be defined using the \code{BinaryMatrix} constructor. In
  the simplest case, you simply pass in the binary data matrix that you
  want to visualize, and the \code{BinaryMatrix is constructed using the matrix's
  existing column and row names.}
}
\usage{
BinaryMatrix(binmat, columnInfo, rowInfo)
}
\arguments{
  \item{binmat}{A binary \code{matrix} of class \code{numeric} or \code{integer}
      with values 0 and 1.}
  \item{columnInfo}{A \code{data.frame} of at least one column containing column
      names. If no \code{colnames} are specified, the column names associated with
      the \code{data} matrix will be used.}
  \item{rowInfo}{A \code{dataframe} of at least one column containing row
      names. If no \code{rownames} are specified, the row names associated with
      the \code{data} matrix will be used.}
}
\section{Slots}{
  \describe{
      \item{\code{@binmat}:}{Object of class \code{matrix}; the binary data 
          used for visualization.}
      \item{\code{@columnInfo}:}{Object of class \code{data.frame}; names and
          definitions of columns.}
      \item{\code{@rowInfo}:}{Object of class \code{data.frame}; names and
          definitions of rows.}
      \item{\code{@info}:}{Object of class \code{list}; identifies \code{$notUsed}
          and \code{$redundant} features.}
      \item{\code{@history}:}{Object of class \code{"character"}; returns a 
          history of manipulations by \code{Mercator} functions to the \code{BinaryMatrix} 
          object, including "Newly created," "Subsetted," "Transposed," "Duplicate 
          features removed," and "Threshed."}
  }
}
\section{Methods}{
  \describe{
      \item{\code{[]}:}{Subsetting by \code{[]} returns a subsetted binary matrix, 
          including subsetted row and column names. Calling \code{@history} will 
          return the history "Subsetted."}
      \item{\code{\link{dim}}:}{returns the dimensions of the \code{@binmat} component 
      of the \code{binaryMatrix} object.}
      \item{\code{\link{summary}}:}{For a given \code{BinaryMatrix}, returns object class,
          dimensions of the \code{@binmat} component, and \code{@history}.}
      \item{\code{\link{t}}:}{Transposes the \code{@binmat} and its associated \code{rowInfo}
          and \code{columnInfo}. Calling \code{@history} will return the history 
          "Subsetted."}
  }
}
\value{
The \code{BinaryMatrix} function returns a new object of \code{binaryMatrix} class.
}
\author{Kevin R. Coombes <krc@silicovore.com}
\note{
Attempting to construct or manipulate a \code{BinaryMatrix} containing NAs, missing values, 
or columns containing exclusively 0 values may introduce error.
}
\seealso{The \code{\link{removeDuplicateFeatures}} function can be used to remove
duplicate columns from the binaryMatrix class before threshing or visualization. 
The \code{\link{threshLGF}} can be used to identify and remove uninformative features
before visualization or further analysis.}
\examples{
my.matrix <- matrix(rbinom(50*100, 1, 0.15), ncol=50)
my.rows <- as.data.frame(paste("R", 1:100, sep=""))
my.cols <- as.data.frame(paste("C", 1:50, sep=""))
my.binmat <- BinaryMatrix(my.matrix, my.cols, my.rows)
summary(my.binmat)
my.binmat <- my.binmat[1:50, 1:30]
my.binmat <- t(my.binmat)
dim(my.binmat)
my.binmat@history
}
\keyword{ cluster }