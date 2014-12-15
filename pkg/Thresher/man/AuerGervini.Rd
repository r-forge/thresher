\name{AuerGervini}
\alias{AuerGervini}
\alias{agDimension}
\title{
Estimating Number of Principal Components Usingthe Auer-Gervini Method
}
\description{
Auer and Gervini [1] dsescribed a Bayesian graphical method for
estimating the number of statisticailly significant principal
components.  We have implemented their method and enhanced it by
automating the final selection.
}
\usage{
AuerGervini(Lambda, epsilon = 2e-16)
agDimension(object)
}
\arguments{
  \item{Lambda}{}
  \item{epsilon}{}
  \item{object}{}
}
\details{
stuff
}
\value{
  stuff
}
\references{

}
\author{Kevin R. Coombes <krc@silicovore.com>}
\seealso{
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Lambda, epsilon = 2e-16) 
{
    if (inherits(Lambda, "SamplePCA")) {
        Lambda <- Lambda@variances
    }
    if (epsilon > 0) {
        Lambda <- Lambda[Lambda > epsilon]
    }
    Lambda <- rev(sort(Lambda))
    rg <- (1:length(Lambda)) - 1
    lowerBounds <- sapply(rg, thetaLower, Lambda = Lambda)
    upperBounds <- sapply(rg, thetaUpper, Lambda = Lambda)
    dLevels <- rev(rg[lowerBounds <= upperBounds])
    changePoints <- rev(lowerBounds[lowerBounds <= upperBounds])[-1]
    new("AuerGervini", Lambda = Lambda, dLevels = dLevels, changePoints = changePoints, 
        lowerBounds = lowerBounds, upperBounds = upperBounds)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
