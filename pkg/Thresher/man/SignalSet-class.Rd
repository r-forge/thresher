\name{SignalSet-class}
\Rdversion{1.1}
\docType{class}
\alias{SignalSet-class}
\alias{SignalSet}

\title{Class \code{"SignalSet"}}
\description{
  We use the term "(continuous) signal" to refer to a weighted sum (by
  default, the mean) of gene-features.  By dichotomizing a continuous
  signals, we obtain a "binary signal". The \code{SignalSet} class
  represents the set of continuous and binary signals obtained after
  clustering the features in a data set.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SignalSet",
  ...)}.  However, users are styrongly discouraged from contructing a
\code{SignalSet} manually. They are only used in the code internal to the
construction of a \code{Reaper} object.
}
\section{Slots}{
  \describe{
    \item{\code{members}:}{Object of class \code{"list"}. Each member of
      the list is a character vector enumerating the features defining
      each signal.}
    \item{\code{continuous}:}{A matrix where the number of columns
      equals the length of the \code{members} list; each column contains
      the mean expression of the (assumed standardized) corresponding features.}
    \item{\code{binary}:}{A matrix where the number of columns
      equals the length of the \code{members} list; each column contains
      expression values dichotmoized to 0 or 1 by splitting the
      conmtinuous siognal at zero.}
    \item{\code{continuousClusters}:}{Object of class \code{"hclust"}
      obtained by clustering samples based on the continuous signals.}
    \item{\code{binaryClusters}:}{Object of class \code{"hclust"}
      obtained by clustering samples based on the binary signals.}
  }
}
\section{Methods}{
No methods defined with class "SignalSet" in the signature.
}
\author{
Kevin R. Coombes <krc@silicovore.com>
}
\note{
  The length of \code{members} and thus the number of signals may be
  smaller than expected from the number of clusters found by
  \code{Reaper}. The implementation of the \code{SignalSet} tries to
  determine if two signals are pointing in opposite directions, which
  could happen if they are postively and negatively correlated
  sets. This behavior is likely to change in the future.
}
\seealso{
  \code{\link{Reaper}}
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
# now extract the signal set
ss <- reap@signalSet
dim(ss@continuous)
dim(ss@binary)
table(ss@binary[,1], ss@binary[,2])
plot(ss@continuousClusters)
}
\keyword{classes}
