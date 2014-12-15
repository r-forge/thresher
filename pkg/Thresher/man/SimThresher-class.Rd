\name{SimThresher-class}
\Rdversion{1.1}
\docType{class}
\alias{SimThresher-class}
\alias{image,SimThresher-method}
\alias{makeFigures,SimThresher-method}

\title{Class \code{"SimThresher"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SimThresher", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{nSample}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{covariance}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{rho}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{data}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{spca}:}{Object of class \code{"SamplePCA"} ~~ }
    \item{\code{loadings}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{gc}:}{Object of class \code{"hclust"} ~~ }
    \item{\code{pcdim}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{delta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{ag}:}{Object of class \code{"AuerGervini"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{Thresher}"}, directly.
}
\section{Methods}{
  \describe{
    \item{image}{\code{signature(x = "SimThresher")}: ... }
    \item{makeFigures}{\code{signature(object = "SimThresher")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("SimThresher")
}
\keyword{classes}
