\name{SimThresher}
\alias{SimThresher}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simt
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SimThresher(ss, nSample, nm = deparse(substitute(ss)), rho = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ss}{
%%     ~~Describe \code{ss} here~~
}
  \item{nSample}{
%%     ~~Describe \code{nSample} here~~
}
  \item{nm}{
%%     ~~Describe \code{nm} here~~
}
  \item{rho}{
%%     ~~Describe \code{rho} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ss, nSample, nm = deparse(substitute(ss)), rho = NULL, 
    ...) 
{
    if (is.null(rho)) {
        rho <- sort(unique(abs(ss[upper.tri(ss)])))[-1]
    }
    require(MASS)
    nFeature <- ncol(ss)
    mu <- rep(0, nFeature)
    simdata <- mvrnorm(nSample, mu, ss)
    colnames(simdata) <- paste("Pr", 1:ncol(simdata), sep = "")
    new("SimThresher", Thresher(simdata, nm, ...), nSample = nSample, 
        covariance = ss, rho = rho)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
