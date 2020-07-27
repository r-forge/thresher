\name{Image}
\alias{RcytoGPS}
\alias{Idiogram Graphs}
\title{Idiogram Graphs}
\description{
  The purpose of this function is to have a visual of any cytogentic event from any dataset as long as you specify the column name, color for the bar graph, and the chromosome number you would want to analyze. To run this function and the other idiogram functions, you have to have data that is formated for Idiograms. I.E. a dataframe that has the summed columns of the LGF mode, a column that specifies what chromosome, along with a column for the chromosome start location, a column for the end location, the arm of the chromosome that the event is taking place on, and the cytogenetic band. The dataframe also has to be transformed into an S4 object using the function CytobandData.
}
\usage{
image(bandData, what = "Loss", chr = 3, pal = "blue", horiz = FALSE)
}
\arguments{
  \describe{
    \item{\code{bandData}}{A dataframe that has the columns of Chromosome,  the start location of the chromosome, the end location of the chromosome, band of the chromosome,  the chromosomal arm, and the column sums of Loss, Gain, Fusion}
    \item{\code{What}}{The column(s) from the data frame you want to analyze. What has to be the column name from the data set. I.E. "A.Loss", "B.Loss", "Loss", "Gain". What can also accept several data columns in a form of a list to do data comparisons on a chromosome.}
    \item{\code{chr}}{ Chr is the specific chromosome you want to analyze, if you want to look at all the chromosomes at once then you can use the command chr = "all"}
    \item{\code{pal}}{ Pal is the arguement if you want to change the color of the idiogram grams.}
    \item{\code{horiz}}{ horiz is the parameter that rotates graph along the x and y axis. By default horiz is FALSE so the images image generates is vertical, but if horiz is TRUE then the graph is rotated to be on the x axis instead of the y axis.}
  }
}
\value{
    Returns a nested list of the binary LGF model as well as a list of which karyotype had a 'success' or 'fail' cytogenetic event. In an event of a 'fail', it could suggestion how to fix the syntax of the inputted karyotype.

}
\source{

}
\references{

}
\examples{
 y <- readLGF("results.JSON", wd)
 x <- readLGF(wd)
 cytoData <- data.frame(x[["CL"]], x[["frequency"]])
 bandData <- CytobandData(cytoData)
 datacolumns <- names(x[["frequency"]])
 image(bandData, what = datacolumns[1:3], chr = 2)

}

\keyword{datasets}



