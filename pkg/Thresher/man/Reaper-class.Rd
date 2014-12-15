\name{Reaper-class}
\Rdversion{1.1}
\docType{class}
\alias{Reaper-class}
\alias{getColors,Reaper-method}
\alias{getSplit,Reaper-method}
\alias{makeFigures,Reaper-method}

\title{Class \code{"Reaper"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Reaper", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{useLoadings}:}{Object of class \code{"logical"} ~~ }
    \item{\code{keep}:}{Object of class \code{"logical"} ~~ }
    \item{\code{nGroups}:}{Object of class \code{"number or miss"} ~~ }
    \item{\code{fit}:}{Object of class \code{"fit or miss"} ~~ }
    \item{\code{allfits}:}{Object of class \code{"list"} ~~ }
    \item{\code{bic}:}{Object of class \code{"number or miss"} ~~ }
    \item{\code{metric}:}{Object of class \code{"character"} ~~ }
    \item{\code{signalSet}:}{Object of class \code{"SignalSet"} ~~ }
    \item{\code{maxSampleGroups}:}{Object of class \code{"numeric"} ~~ }
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
    \item{getColors}{\code{signature(object = "Reaper")}: ... }
    \item{getSplit}{\code{signature(object = "Reaper")}: ... }
    \item{makeFigures}{\code{signature(object = "Reaper")}: ... }
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
showClass("Reaper")
}
\keyword{classes}
