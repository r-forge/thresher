setClass("SimThresher",
         representation=list("Thresher",
           nSample="numeric",
           covariance="matrix",
           rho="numeric"
           ))

setMethod("image", "SimThresher", function(x, ...) {
  image(x@covariance, ...)
})

SimThresher <- function(ss, nSample, nm=deparse(substitute(ss)), rho=NULL,
                        agfun=agDimTwiceMean, ...) {
  if (is.null(rho)) {
    rho <- sort(unique(abs(ss[upper.tri(ss)])))[-1]
  }
  nFeature <- ncol(ss)
  mu <- rep(0, nFeature)
  simdata <- mvrnorm(nSample, mu, ss)
  colnames(simdata) <- paste("Pr", 1:ncol(simdata), sep='')
  new("SimThresher", Thresher(simdata, nm, agfun=agfun, ...),
      nSample=nSample,
      covariance = ss,
      rho=rho)
}

setMethod("makeFigures", "SimThresher", function(object, DIR=NULL, ...) {
  if (is.null(DIR)) opar <- par(ask=TRUE)
  fname <- gsub("\\.", "-", object@name) # latex-safe
# fig0
  if (!is.null(DIR)) {
    png(filename=file.path(DIR, paste(fname, "-00-ss.png", sep="")),
        width=640, height=640, bg="white")
  }
  image(object, col=blueyellow(64), zlim=c(-1,1))
  if(!is.null(DIR)) dev.off()
  callNextMethod()
  invisible(object)
})

