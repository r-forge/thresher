\name{generics}
\alias{makeFigures}
\alias{makeFigures,ANY-method}
\alias{scatter}
\alias{scatter,ANY-method}
\alias{heat}
\alias{heat,ANY-method}
\alias{getSplit}
\alias{getSplit,ANY-method}
\alias{getColors}
\alias{getColors,ANY-method}
\alias{getStyles}
\alias{getStyles,ANY-method}
\title{Generic Methods in the Thresher package}
\description{
  New generic functions for threshing and reaping datasets.
}
\usage{
\S4method{getColors}{ANY}(object, \dots)
\S4method{getSplit}{ANY}(object, \dots)
\S4method{getStyles}{ANY}(object, \dots)
\S4method{scatter}{ANY}(object, \dots)
\S4method{heat}{ANY}(object, \dots)
\S4method{makeFigures}{ANY}(object, DIR=NULL, \dots)
}
\arguments{
  \item{object}{an object of any class, but typically
    \code{\link{Thresher}}, \code{\link{Reaper}}, or
    \code{\link{SimThresher}}. } 
  \item{DIR}{a character string containing the name of an existing directory.}
  \item{\dots}{extra arguments for generic routines.}
}
\details{
  The methods \code{getColors}, \code{getSplit}. and \code{getStyles}
  are used to extract the implicit color assignments used in standard
  plots of objects of the \code{Thresher} or \code{Reaper} classes.

  The \code{heat} and \code{scatter} methods represent generic heat maps
  and scatter plots, respectively.

  The \code{makeFigures} method is to generate a standard suite of
  figures for an object.  If the \code{DIR} argument is not NULL, then
  the figures will be written to the indicated directory in PNG
  format. Otherwise, the figures will be displayed interactivey on
  screen, waiting for user input to show each plot.
}
\value{
  The form of the value returned by these functions may change depending
  on the class of its argument. See the documentation of the particular
  methods for details of what is produced by that method. 

  The \code{getColors} method should return a vector of colors with
  length equal to the number of columns in a data set.

  The \code{getSplit} method should return a vector of colors with
  length equal to the number of rows in a data set.

  The \code{heat}, \code{scatter}, and \code{makeFigures} methods are
  called for their side effects of producing plots.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com}
}
\seealso{
  \code{\link{Thresher}}, \code{\link{Reaper}}, and \code{\link{SimThresher}}.
}
\keyword{methods}

