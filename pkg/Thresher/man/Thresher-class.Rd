\name{Thresher-class}
\Rdversion{1.1}
\docType{class}
\alias{Thresher-class}
\alias{getColors,Thresher-method}
\alias{getSplit,Thresher-method}
\alias{heat,Thresher-method}
\alias{makeFigures,Thresher-method}
\alias{plot,Thresher,missing-method}
\alias{scatter,Thresher-method}
\alias{screeplot,Thresher-method}

\title{Class \code{"Thresher"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Thresher", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
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
\section{Methods}{
  \describe{
    \item{getColors}{\code{signature(object = "Thresher")}: ... }
    \item{getSplit}{\code{signature(object = "Thresher")}: ... }
    \item{heat}{\code{signature(object = "Thresher")}: ... }
    \item{makeFigures}{\code{signature(object = "Thresher")}: ... }
    \item{plot}{\code{signature(x = "Thresher", y = "missing")}: ... }
    \item{scatter}{\code{signature(object = "Thresher")}: ... }
    \item{screeplot}{\code{signature(x = "Thresher")}: ... }
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
showClass("Thresher")
}
\keyword{classes}
