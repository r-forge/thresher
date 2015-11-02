#######################################################
library(PCDimension)
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
# test external interface on random data

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

rm(spca)
######################################################
# get the previously computed simulated SamplePCA object
data(spca)

vars <- spca@variances
obj <- AuerGervini(vars, dd=dim(spca@scores))
summary(obj)

plot(obj)

f <- makeAgCpmFun("Exponential")
plot(obj, list(agDimTwiceMean, agDimKmeans,
               agDimTtest, agDimTtest2, agDimCPT, f))

