\name{CytobandData-methods}
\alias{image,CytobandData-method}
\alias{barplot,CytobandData-method}
\alias{Idiogram Graphs}
\title{Idiogram Graphs}
\description{
  The purpose of this function is to have a visual of any cytogentic event from any dataset as long as you specify the column name, color for the bar graph, and the chromosome number you would want to analyze. To run this function and the other idiogram functions, you have to have data that is formated for Idiograms. I.E. a dataframe that has the summed columns of the LGF mode, a column that specifies what chromosome, along with a column for the chromosome start location, a column for the end location, the arm of the chromosome that the event is taking place on, and the cytogenetic band. The dataframe also has to be transformed into an S4 object using the function CytobandData.
}
\usage{
\S4method{image}{CytobandData}(x, chr, what,
         pal = palette(), nrows = 2, labels = NULL,
         horiz = FALSE, axes = chr != "all", debug = FALSE, legend = FALSE,
         sigcolumn = NA, sigcut = 0.01, alpha = 63)
\S4method{barplot}{CytobandData}(height, what, col = "blue", altcol = "#FED4C4",
                   ylab = "Percent", hline = NULL,
                   xform = function(x) x, ...)
}
\arguments{
  \item{x}{An pbject of the \code{\link{CytobandData-class}}, which
    combines numerical data at cytoband resolution with information on
    the chromosomal locations of the cytobands.}
  \item{chr}{The specific chromosome you want to view, typically in
    \code{c(1:22, "X", "Y")}. if you want to see all the chromosomes
    at once then you can set \code{chr = "all"}.}
  \item{what}{A vector or list containing the names of the numerical
    column(s) that you want to display from the data frame. The
    \code{barplot} method only shows a single data column at a time. The
    plots resulting from the \code{image} method change depending on
    whether you supply a vector or a list, as well as on the length. For
    more details, see the vignettes, espefcialy the image gallery.}
  \item{pal}{a charcter vector containing the colors you want to use
    for different data shown in the plot.}
  \item{horiz}{A logical value determining whether images present the
    chromosome idiograms horizontally or vertically.}
  \item{nrows}{Only used when \code{chr = "all"} to determine the number
    of rows to use to display different chromosomes. Must be an integer
    between 1 and 4.}
  \item{labels}{Only use when \code{what} is a character vector (not a
    list) and \code{chr} is mot equal to \code{"all"}. Used to label
    differnt columns of displayed data.}
  \item{axes}{Logical value; should axes be displayed?}
  \item{legend}{Logicla value; should a legend be added to the plot.}
  \item{sigcolumn}{The three parameters \code{sigcolumn},
    \code{sigcut}, and \code{alpha} are used as a set. The first names a
    numerical column in the data set used to define "significance". The
    second is a vector of cutoffs that mark levels of significance. The
    latter is a number between 0 and 255 denoting the transparency level
    assigned to the color for each significance level.}
  \item{sigcut}{See \code{sigcolumn}.}
  \item{alpha}{See \code{sigcolumn}.}
  \item{debug}{Logical value; should the method print out debugging
    information. Probably best to ignore.}

  \item{height}{An pbject of the \code{\link{CytobandData-class}}, which
    combines numerical data at cytoband resolution with information on
    the chromosomal locations of the cytobands.}
  \item{col}{a character vector containing the colors you want to use
    for different data shown in the plot.}
  \item{altcol}{Determines the color used, lternating with white, in the
    x-axis plot of all chromosomes.}
  \item{ylab}{Label for the y-axis in a barplot.}
  \item{hline}{Numeric vector of heights at which toadd a horizontal
    line.}
  \item{xform}{Function to transform the data before plotting. Defauly
    is the identity map, which does nothing.}
  \item{\dots}{Ignored.}
}
\value{
    Returns a nested list of the binary LGF model as well as a list of which karyotype had a 'success' or 'fail' cytogenetic event. In an event of a 'fail', it could suggestion how to fix the syntax of the inputted karyotype.

}
\references{
  Abrams ZB, Tally DG, Coombes KR. RCytoGPS: An R Package for Analyzing
  and Visualizing Cytogenetic Data. In preparation.

  Abrams ZB, Tally DG,  Zhang L, Coombes CE, Payne PRO, Abruzzo LV,
  Coombes KR. Pattern recognition in lymphoid malignancies using CytoGPS
  and Mercator. Under review.
}
\examples{
jsonDir <- system.file("Examples/JSONfiles", package = "RCytoGPS")
x <- readLGF(folder = jsonDir)
cytoData <- data.frame(x[["CL"]], x[["frequency"]])
bandData <- CytobandData(cytoData)
datacolumns <- names(x[["frequency"]])
barplot(bandData, what = datacolumns[1], col="forestgreen")

image(bandData, what = datacolumns[1:3], chr = 2)
image(bandData, what = datacolumns[1:3], chr = "all")
image(bandData, what = as.list(datacolumns[1:2]), chr = 2)
image(bandData, what = as.list(datacolumns[1:2]), chr = "all")
}
\keyword{graphics}



