\name{threshLGF}
\alias{threshLGF}
\alias{ThreshedBinaryMatrix}
\alias{ThreshedBinaryMatrix-class}
\title{
Threshing and Reaping the \code{BinaryMatrix}
}
\description{
The \code{\link{threshLGF}} function produces an object of class 
\code{\link{ThreshedBinaryMatrix}} from threshing on an object of class 
\code{\link{BinaryMatrix}}. 

The function \code{\link{threshLGF}} and the \code{\link{ThreshedBinaryMatrix}} 
object can be used to access the functionality of the \code{\link{Thresher}} 
R-package within Mercator.
}
\usage{
threshLGF(object, cutoff = 0)
}
\arguments{
  \item{object}{An object of class {BinaryMatrix}
  }
  \item{cutoff}{The value of \code{delta} set to demarcate an 
  uninformative feature. Features with a value greater than the 
  cutoff will be kept.
  }
}
\details{
  The \code{\link{Thresher}} R-package provides a variety of functionalities
  for data filtering and the identification of and reduction to "informative" features.
  It performs clustering using a combination of outlier detection, principal
  component analysis, and von Mises Fisher mixture models. By identifying 
  significant features, Thresher performs feature reduction through the 
  identification and removal of noninformative features and the nonbiased 
  calculation of the number of groups (K) for down-stream use.
}
\value{
  \code{threshLGF} returns an object of class \code{ThreshedBinaryMatrix}.
  The \code{ThreshedBinaryMatrix} object retains all the functionality, 
  slots, and methods of the \code{\link{BinaryMatrix}} object class with added 
  features. After threshing, the \code{ThreshedBinaryMatrix} records the 
  \code{history}, "Threshed."
}
\section{Additional Slots}{
  \describe{
    \item{\code{thresher}:  }{Returns the functions of the \code{Thresher} 
    object class of the \code{\link{Thresher}} R-package.
    }
    \item{\code{reaper}: }{Returns the functions of the \code{Reaper} 
    object class of the \code{\link{Thresher}} R-package.
    }
  }
}
\references{
Wang, M., Abrams, Z. B., Kornblau, S. M., & Coombes, K. R. (2018). Thresher: 
determining the number of clusters while removing outliers. BMC bioinformatics, 
19(1), 9.
}
\author{Kevin R. Coombes <krc@silicovore.com}
\note{
The \code{\link{Thresher}} R-package applies the Auer-Gervini statistic 
for principal component analysis, outlier detection, and identification 
of uninformative features on a \code{matrix} of class \code{integer} or 
\code{numeric}.

An initial \code{delta} of 0.3 is recommended.
}
\seealso{
The \code{threshLGF} function creates a new object of class 
\code{ThreshedBinaryMatrix} from an object of class \code{BinaryMatrix}.
}
\examples{
#Create a BinaryMatrix
set.seed(52134)
my.matrix <- matrix(rbinom(50*100, 1, 0.15), ncol=50)
my.rows <- as.data.frame(paste("R", 1:100, sep=""))
my.cols <- as.data.frame(paste("C", 1:50, sep=""))
my.binmat <- BinaryMatrix(my.matrix, my.cols, my.rows)
summary(my.binmat)

#Identify delta cutoff and thresh
my.binmat <- threshLGF(my.binmat)
Delta <- my.binmat@thresher@delta
sort(Delta)
hist(Delta, breaks=15, main="", xlab="Weight")
abline(v=0.3, col='red')
my.binmat <- threshLGF(my.binmat, cutoff = 0.3)
summary(my.binmat)

#Principal Component Analysis
my.binmat@reaper@pcdim
my.binmat@reaper@nGroups
#namespace issue: plot(my.binmat@reaper@ag)
#abline(h=1, col="red")
#screeplot(my.binmat@reaper)
#abline(v=6, col="green")
  }
}
\keyword{ cluster }
