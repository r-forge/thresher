\name{matchLabels}
\alias{matchLabels}
\alias{labelMatcher}
\alias{labelAccuracy}
\alias{countAgreement}
\alias{bestMetric}
\title{
Match Arbitrary Class Assignments Across Methods
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
matchLabels(tab)
labelMatcher(tab, verbose = FALSE)
labelAccuracy(data, labels)
countAgreement(tab)
bestMetric(data, labels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{
%%     ~~Describe \code{tab} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
factor1 <- sample(c("A", "B", "C"), 30, replace=TRUE)
factor2 <- rep(c("X", "Y", "Z"), each=10)
tab <- table(factor1, factor2)
matchLabels(tab)
labelMatcher(tab)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ classif }
\keyword{ multivariate }
