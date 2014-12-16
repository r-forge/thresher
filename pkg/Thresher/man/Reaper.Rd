\name{Reaper}
\alias{Reaper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reaper
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Reaper(thresher, useLoadings = FALSE, cutoff = 0.3,
       metric = NULL, linkage="ward.D2",
       maxSampleGroups = 0, verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{thresher}{
%%     ~~Describe \code{thresher} here~~
}
  \item{useLoadings}{
%%     ~~Describe \code{useLoadings} here~~
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
}
  \item{metric}{
%%     ~~Describe \code{metric} here~~
}
  \item{LINKAGE}{
%%     ~~Describe \code{metric} here~~
}
  \item{maxSampleGroups}{
%%     ~~Describe \code{maxSampleGroups} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (thresher, useLoadings = FALSE, cutoff = 0.3, metric = NULL, 
    maxSampleGroups = 0, verbose = TRUE, ...) 
{
    if (verbose) 
        cat(thresher@name, "\n", file = stderr())
    keep <- thresher@delta > cutoff
    m <- ifelse(is.null(metric), "pearson", metric)
    cleaned <- Thresher(thresher@data[, keep], paste(thresher@name, 
        "cleaned", sep = "."), metric = m, ...)
    tab <- 0
    counter <- 0
    while (any(tab == 0) & counter < 5) {
        counter <- counter + 1
        fits <- .fitModels(cleaned, useLoadings)
        if (length(fits) == 0) 
            next
        bic <- sapply(fits, BIC)
        woo <- which(bic == min(bic))
        ng <- as.integer(sub("NC=", "", names(woo)))
        fit <- fits[[woo]]
        gassign <- factor(predict(fit), levels = 1:ng)
        tab <- table(gassign)
    }
    if (length(fits) == 0) {
        bic <- ng <- fit <- NA
        metric <- "no fit"
        sigset <- new("SignalSet")
    }
    else {
        if (is.null(metric)) {
            pp <- factor(paste("G", predict(fit), sep = ""))
            metric <- bestMetric(cleaned@data, pp)
            cleaned@gc <- hclust(distanceMatrix(cleaned@data, 
                metric, p = 1), "ward")
        }
        if (any(tab == 0)) {
            sigset <- new("SignalSet")
        }
        else {
            sigset <- .findSignals(cleaned, fit, ng)
        }
    }
    new("Reaper", cleaned, useLoadings = useLoadings, keep = keep, 
        nGroups = ng, fit = fit, allfits = fits, bic = bic, metric = metric, 
        signalSet = sigset, maxSampleGroups = maxSampleGroups)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
