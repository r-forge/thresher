\name{rnames}
\alias{RcytoGPS}
\title{rnames}
\description{
   The purpose of this function is to recursively extract the element index of karyotypes from a list using a specific keyword found within the list to generate an ID for each karyotype. The function can also produce additional columns such as the clones which allows the function readLGF to check if there is a disagreement in the clones, a duplicate in the elements, and the length of each element.
}
\usage{
rnames(list, item)
}
\arguments{
  \describe{
    \item{\code{list}}{A list that contains the loss gain fusion binary model of each karyotype and their clones.list could also be in a form of a nested list.}
    \item{\code{item}}{Item is a specific keyword that has to be able to be found in the given list. item is what rnames look for in the list to extract the element index.}
  }
}
\value{
  Returns an array of integers where a single elements is expressed in a set of 4. You can change the array into a data frame and change the dimension to have a more visually appealing format.
}
\source{

}
\references{

}
\examples{
 rnames(lgf, "loss")
 x <- as.data.frame(matrix(rnames(lgf, "loss"), ncol = 4, byrow = TRUE))


}
\keyword{datasets}

