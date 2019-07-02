\name{Mercator}
\alias{Mercator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Mercator(binaryMat, metric, method, K, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{binaryMat}{
%%     ~~Describe \code{binaryMat} here~~
}
  \item{metric}{
%%     ~~Describe \code{metric} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
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
function (binaryMat, metric, method, K, ...) 
{
    DistMat <- binaryDistance(binaryMat@binmat, metric)
    M <- attr(DistMat, "comment")
    if (!is.null(M)) 
        metric <- M
    poof <- pam(DistMat, k = K, diss = TRUE, cluster.only = TRUE)
    dispSet <- makeDisplay(poof)
    colv <- dispSet$colv
    symv <- dispSet$symv
    names(colv) <- names(symv) <- labels(DistMat)
    view <- list()
    ob <- new("Mercator", metric = metric, distance = DistMat, 
        view = view, colv = colv, symv = symv)
    addVisualization(ob, method, ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
