\name{Thresher}
\alias{Thresher}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Thresher
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Thresher(data, nm = deparse(substitute(data)), FUZZ = 0.005,
         metric = "pearson", lonkage="ward.D2",
         method = c("broken.stick", "auer.gervini"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{nm}{
%%     ~~Describe \code{nm} here~~
}
  \item{FUZZ}{
%%     ~~Describe \code{FUZZ} here~~
}
  \item{metric}{
%%     ~~Describe \code{metric} here~~
}
  \item{linkage}{
%%     ~~Describe \code{metric} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (data, nm = deparse(substitute(data)), FUZZ = 0.005, 
    metric = "pearson", method = c("broken.stick", "auer.gervini")) 
{
    std <- scale(data)
    spca <- SamplePCA(t(std))
    ag <- AuerGervini(spca)
    method <- match.arg(method)
    pcdim <- switch(method, broken.stick = bsDimension(spca), 
        auer.gervini = agDimension(ag))
    deltaDim <- max(1, pcdim)
    lambda <- sqrt(spca@variances)
    loadings <- sweep(spca@components, 2, lambda, "*")
    delta <- sqrt(apply(loadings[, 1:deltaDim, drop = FALSE]^2, 
        1, sum))
    gc <- hclust(distanceMatrix(std, metric), "ward")
    new("Thresher", name = nm, data = data, spca = spca, loadings = loadings, 
        gc = gc, pcdim = pcdim, delta = delta, ag = ag)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
