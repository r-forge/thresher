#######################################################
# INTERNAL
# functions based on descriptions in Auer-Gervini, 2008

fhat <- function(d, Lambda) {
  m <- length(Lambda)
  use <- Lambda[(d+1):m]
  geo <- exp(mean(log(use)))
  ari <- mean(use)
  (m-d)*log(geo/ari)
}

dhat <- function(theta, Lambda) {
  ds <- seq(0, length(Lambda)-1)
  vals <- sapply(ds, function(d, theta, Lambda) {
    fhat(d, Lambda) + theta*(length(Lambda)-1-d)
  }, theta=theta, Lambda=Lambda)
  max(ds[which(vals==max(vals))])
}

thetaLower <- function(d, Lambda) {
  m <- length(Lambda)
  if(d == m - 1) return(0)
  fd <- fhat(d, Lambda)
  kset <- (d+1):(m-1)
  temp <- sapply(kset, function(k) {
    (fhat(k, Lambda) - fd) / (k - d)
  })
  max(temp)
}

thetaUpper <- function(d, Lambda) {
  if (d == 0) return(Inf)
  m <- length(Lambda)
  fd <- fhat(d, Lambda)
  kset <- 0:(d-1)
  temp <- sapply(kset, function(k) {
    (fhat(k, Lambda) - fd) / (k - d)
  })
  min(temp)
}

#######################################################
# EXTERNAL

# broken-stick model
brokenStick <- function(k, n) {
    if (any(n < k)) stop('bad values')
    x <- 1/(1:n)
    sapply(k, function(k0) sum(x[k0:n])/n)
}

# broken stick to get dimension
bsDimension <- function(lambda, FUZZ=0.005) {
  if(inherits(lambda, "SamplePCA")) {
    lambda <- lambda@variances
  }
  N <- length(lambda)
  bs  <- brokenStick(1:N, N)
  fracVar <- lambda/sum(lambda)
  which(fracVar - bs < FUZZ)[1] - 1
}

#######################################################
# S4 interface

setClass("AuerGervini",
         representation=list(
           Lambda="numeric",
           dimensions="numeric",
           dLevels="numeric",
           changePoints="numeric",
           lowerBounds="numeric",
           upperBounds="numeric"
           ))

AuerGervini <- function(Lambda, dd=NULL, epsilon=2e-16) {
  if (inherits(Lambda, "SamplePCA")) {
    dd <- dim(Lambda@scores)
    Lambda <- Lambda@variances
  }
  if (epsilon > 0) {
    Lambda <- Lambda[Lambda > epsilon]
  }
  Lambda <- rev(sort(Lambda))
  rg <- (1:length(Lambda)) - 1
  lowerBounds <- sapply(rg, thetaLower, Lambda=Lambda)
  upperBounds <- sapply(rg, thetaUpper, Lambda=Lambda)
  dLevels <- rev(rg[lowerBounds <= upperBounds])
  changePoints <- rev(lowerBounds[lowerBounds <= upperBounds])[-1]
  new("AuerGervini", Lambda=Lambda, dimensions = dd,
      dLevels=dLevels, changePoints=changePoints,
      lowerBounds=lowerBounds, upperBounds=upperBounds)
}

estimateTop <- function(object) {
  oldAndCrude <- -2*log(0.01)/length(object@Lambda) # in case we want to revert
  n <- object@dimensions[1] # nrows
  m <- object@dimensions[1] # ncolumns
  delta <- 1 # why do we do it this way?
  magic <- 18.8402+1.9523*m-0.0005*m^2 # from a linear model fit on simulated data
  modelBased <- ifelse(n >= m,
                       magic/n,
                       magic*n/m^2) # please explain why
  max(oldAndCrude,
      modelBased,
      1.03*object@changePoints)
}

agDimension <- function(object) {
  stepLength <- diff(c(object@changePoints, estimateTop(object)))
  if (length(stepLength) > 3) {
    magic <- (stepLength > 2*mean(stepLength))
  } else {
    magic <- (stepLength == max(stepLength))
  }
  object@dLevels[1+which(magic)[1]]
}

setMethod("plot", c("AuerGervini", "missing"), function(x, y,
        main="Bayesian Sensitivity Analysis",
        ...) {
  top <-estimateTop(x)
  fun <- stepfun(x@changePoints, x@dLevels)
  plot(fun, xlab="Prior Theta", ylab="Number of Components",
     main=main, xlim=c(0, top), ...)
 abline(h=agDimension(x), lty=2, lwd=2, col='pink')
 invisible(x)
})

setMethod("summary", "AuerGervini", function(object, ...) {
  cat("An '", class(object), "' object that estimates the number of ",
      "principal components to be ", agDimension(object), ".\n", sep="")
})
