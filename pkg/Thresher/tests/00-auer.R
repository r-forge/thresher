#######################################################
library(Thresher)
# simulate one unstructured dataset
set.seed(992500)
NC <- 15
NS <- 200
ranData <- matrix(rnorm(NS*NC, 6), ncol=NC)
# perform PCA and get the set of variances/eigenvectors
spca <- SamplePCA(ranData)
lam <- spca@variances[1:(NC-1)] # remove 0 eigenvalue
rg <- 0:(NC-2)

######################################################
# Test some internal routines

# Theory says that fhat is non-decreasing
fh <- sapply(rg, Thresher:::fhat, Lambda=lam)
all( diff(fh) >= 0 ) # should be 'TRUE'

plot( rg, fh, type='b', pch=16 )

# Theory says that dhat is non-increasing
thetas <- seq(0, 0.2, by=0.01)
dset <- sapply(thetas, Thresher:::dhat, Lambda=lam)
all( diff(dset) <= 0 ) # should be 'TRUE'

plot(thetas, dset, type='l')

######################################################
# test external interface

# Auer-Gervini Bayesian method
obj <- AuerGervini(spca)
summary(obj)
agDimension(obj)

plot(obj)

# broken-stick model
bs <- brokenStick(1:NC, NC)
bs0 <- brokenStick(1:(NC-1), (NC-1))

pts <- screeplot(spca, ylim=c(0, 0.2))
lines(pts, bs, type='b', col='blue', lwd=2, pch=16)
lines(pts[-NC], bs0, type='b', col='red', lwd=2, pch=16)

# get the previously computed simulated data sets
data(savedSims)

idx <- sample(length(savedSims), 1)
ng <- c(1, 2, 1, 2, 2)[idx]

s1 <- savedSims[[idx]]
vars <- s1@spca@variances
obj <- AuerGervini(vars)
summary(obj)

plot(obj, main=paste(idx, ";  True Groups =", ng))
abline(h=ng, lty=3)

#par(ask=T)
for (i in 1:5) {
  s1 <- savedSims[[i]]
  ng <- c(1, 2, 1, 2, 2)[i]
  obj <- AuerGervini(s1@spca)
  plot(obj)
  abline(h=ng, lty=3)
  nSam <- nrow(s1@data)
  nProtein <- ncol(s1@data)
  simData <- matrix(rnorm(nSam*nProtein), ncol=nProtein)
  spca <- SamplePCA(simData)
  vars <- spca@variances[-nProtein]
  obj <- AuerGervini(vars)
  plot(obj)
  abline(h=0, lty=3, col='pink')
}

