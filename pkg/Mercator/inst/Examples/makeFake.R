library(Umpire)

set.seed(391629)
## Set up survival outcome; baseline is exponential
sm <- SurvivalModel(baseHazard=1/5, accrual=5, followUp=1)
## Build a CancerModel with 8 subtypes
nBlocks <- 20    # number of possible hits
cm <- CancerModel(name="cansim",
                  nPossible = nBlocks,
                  nPattern = 8,
                  OUT = function(n) rnorm(n, 0, 1), 
                  SURV= function(n) rnorm(n, 0, 1),
                  survivalModel=sm)
## Include 100 blocks/pathways that are not hit by cancer
nTotalBlocks <- nBlocks + 30
## Assign values to hyperparameters
## block size
blockSize <- round(rnorm(nTotalBlocks, 15, 4))
## log normal mean hypers
mu0    <- 6
sigma0 <- 1.5
## log normal sigma hypers
rate   <- 28.11
shape  <- 44.25
## block corr
p <- 0.6
w <- 5
## Set up the baseline Engine
rho <- rbeta(nTotalBlocks, p*w, (1-p)*w)
base <- lapply(1:nTotalBlocks,
               function(i) {
                 bs <- blockSize[i]
                 co <- matrix(rho[i], nrow=bs, ncol=bs)
                 diag(co) <- 1
                 mu <- rnorm(bs, mu0, sigma0)
                 sigma <- matrix(1/rgamma(bs, rate=rate, shape=shape), nrow=1)
                 covo <- co *(t(sigma) %*% sigma)
                 MVN(mu, covo)
               })
eng <- Engine(base)
## Alter the means if there is a hit
altered <- alterMean(eng, normalOffset, delta=0, sigma=1)
## Build the CancerEngine using character strings
object <- CancerEngine(cm, "eng", "altered")
## Or build it using the actual Engine components
ob <- CancerEngine(cm, eng, altered)
summary(object)
summary(ob)
## Simulate the data
dset <- rand(object, 300)

fakedata <- dset$data
fakeclin <- dset$clin
colnames(fakeclin)[1] <- "Type"

save(fakedata, fakeclin, file = "../../data/fakedata.rda")

if (FALSE) {
  library(Polychrome)
  colset <- palette36.colors(10)[3:10]
  library(ClassDiscovery)
  spca <- SamplePCA(fakedata)
  plot(spca, split=factor(fakeclin$Type), col = colset, cex=1.4)
  library(Rtsne)
  tis <- Rtsne(t(fakedata))
  plot(tis$Y[, 1], tis$Y[, 2], pch = 16, col=colset[fakeclin$Type], cex=1.4)

  library(Mercator)
  mercury <- Mercator(distanceMatrix(fakedata, "euclid"), "euclid", "hclust", 8)
  mercury <- addVisualization(mercury, "mds")
  mercury <- addVisualization(mercury, "tsne")
  plot(mercury, view = "hclust")
  plot(mercury, view = "mds", cex=1.5)
  plot(mercury, view = "tsne", cex=1.5)
  mercury <- setClusters(mercury, fakeclin$Type)
  plot(mercury, view = "hclust")
  plot(mercury, view = "mds", cex=1.5)
  plot(mercury, view = "tsne", cex=1.5)
}

