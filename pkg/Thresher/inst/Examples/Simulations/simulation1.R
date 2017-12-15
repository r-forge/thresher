
##################################################################################
###               Simulation Study on Computing Number of Clusters             ###
##################################################################################


###------------------------------------------------------------------------------------
### load the libraries

library(PCDimension)
library(nFactors)
library(cluster)
library(Thresher)
library(NbClust)
library(rstiefel)
existing <- file.path("..", "ExistingMethods")
source(file.path(existing, "NbClust.txt")) # remove set.seed(1) lines 
source(file.path(existing, "SCOD_functions.R"))

## set control parameters

N <- 96
p <- 24
K <- 3
gsize <- p/K # size of each block in correlation matrices 5-13


###--------------------------------------------------------------------------------------
### define correlation matrices
### In manuscript, we are using correlation matrices 4 - 19

cormatlist <- list()
for (j in 1:19) {
  cormatlist[[j]] <- matrix(0,p,p)
}

# covariance matrix 1-3
lambdavec <- 1/10*rep(1,p)
for (l in 1:5) {
  lambdavec[l] <- 1/l
}
Tau <- rustiefel(p,p)
cormatlist[[1]] <- cormatlist[[2]] <- cormatlist[[3]] <- Tau%*%diag(lambdavec)%*%t(Tau)

# correlation matrix 4
diag(cormatlist[[4]]) <- 1

# correlation matrix 5
cormatlist[[5]] <- matrix(0.8,p,p)
diag(cormatlist[[5]]) <- 1

# correlation matrix 6
cormatlist[[6]] <- matrix(0.3,p,p)
diag(cormatlist[[6]]) <- 1

# correlation matrix 7
for (k in 1:K) {
  cormatlist[[7]][(k-1)*gsize+1:gsize,(k-1)*gsize+1:gsize] <- 0.8
}
diag(cormatlist[[7]]) <- 1

# correlation matrix 8
for (k in 1:K) {
  cormatlist[[8]][(k-1)*gsize+1:gsize,(k-1)*gsize+1:gsize] <- 0.3
}
diag(cormatlist[[8]]) <- 1

# correlation matrix 9
cormatlist[[9]] <- matrix(0.3,p,p)
for (k in 1:K) {
  cormatlist[[9]][(k-1)*gsize+1:gsize,(k-1)*gsize+1:gsize] <- 0.8
}
diag(cormatlist[[9]]) <- 1

# correlation matrix 10
cormatlist[[10]][1:((K-1)*gsize),1:((K-1)*gsize)] <- 0.1
lam10 <- c(0.8,0.3)
for (k in 1:(K-1)) {
  cormatlist[[10]][(k-1)*gsize+1:gsize,(k-1)*gsize+1:gsize] <- lam10[k]
}
diag(cormatlist[[10]]) <- 1

# correlation matrix 11
lam11 <- c(0.8,0.3)
for (k in 1:(K-1)) {
  cormatlist[[11]][(k-1)*gsize+1:gsize,(k-1)*gsize+1:gsize] <- lam11[k]
}
diag(cormatlist[[11]]) <- 1

# correlation matrix 12
cormatlist[[12]][1:gsize,1:gsize] <- 0.8
diag(cormatlist[[12]]) <- 1

# correlation matrix 13
cormatlist[[13]][1:gsize,1:gsize] <- 0.3
diag(cormatlist[[13]]) <- 1

# correlation matrix 14
cormatlist[[14]] <- matrix(0.8,p,p)
cormatlist[[14]][1:(p/2),p/2+1:(p/2)] <- -0.8
cormatlist[[14]][p/2+1:(p/2),1:(p/2)] <- -0.8
diag(cormatlist[[14]]) <- 1

# correlation matrix 15
cormatlist[[15]] <- matrix(0.3,p,p)
cormatlist[[15]][1:(p/2),p/2+1:(p/2)] <- -0.3
cormatlist[[15]][p/2+1:(p/2),1:(p/2)] <- -0.3
diag(cormatlist[[15]]) <- 1

# correlation matrix 16
cormatlist[[16]][1:(p/2),1:(p/2)] <- 0.8
cormatlist[[16]][1:(p/4),p/4+1:(p/4)] <- -0.8
cormatlist[[16]][p/4+1:(p/4),1:(p/4)] <- -0.8
cormatlist[[16]][p/2+1:(p/2),p/2+1:(p/2)] <- cormatlist[[16]][1:(p/2),1:(p/2)]
diag(cormatlist[[16]]) <- 1

# correlation matrix 17
cormatlist[[17]][1:(p/2),1:(p/2)] <- 0.3
cormatlist[[17]][1:(p/4),p/4+1:(p/4)] <- -0.3
cormatlist[[17]][p/4+1:(p/4),1:(p/4)] <- -0.3

cormatlist[[17]][p/2+1:(p/2),p/2+1:(p/2)] <- cormatlist[[17]][1:(p/2),1:(p/2)]
diag(cormatlist[[17]]) <- 1

# correlation matrix 18
cormatlist[[18]][1:(p/2),1:(p/2)] <- cormatlist[[16]][1:(p/2),1:(p/2)]
cormatlist[[18]][p/2+1:(p/2),p/2+1:(p/2)] <- 0.8
diag(cormatlist[[18]]) <- 1

# correlation matrix 19
cormatlist[[19]][1:(p/2),1:(p/2)] <- cormatlist[[17]][1:(p/2),1:(p/2)]
cormatlist[[19]][p/2+1:(p/2),p/2+1:(p/2)] <- 0.3
diag(cormatlist[[19]]) <- 1


###----------------------------------------------------------------------------------------
### define functions and notations to save results

## define criteria in PCDimension
f <- makeAgCpmFun("Exponential")
agfuns <- list(twice=agDimTwiceMean, specc=agDimSpectral,
               km=agDimKmeans, km3=agDimKmeans3, 
               tt=agDimTtest, tt2=agDimTtest2,
               cpt=agDimCPT, cpm=f)
l <- 7  # firstly use 7th criterion -- CPT
# then use the 1st criterion -- TwiceMean

# options(show.error.messages = FALSE)

M <- 1000  # simulate 1000 datasets
minnc <- 1
maxnc <- 8

## use criterion CPT (first choice) and SCOD
pcdimmat <- matrix(0, M, 19) # save pc components
bicdimmat <- matrix(0, M, 19) # save number of clusters based on BIC
noiselist <- list() # save noise features
scod.ncmat <- matrix(0, M, 19)
scod.noiselist <- list()
for (i in 1:19) {
  noiselist[[i]] <- list()
  scod.noiselist[[i]] <- list()
}

clustlist <- list() # save features' grouping labels
scod.clustlist <- list()
for (i in 1:19) {
  clustlist[[i]] <- list()
  scod.clustlist[[i]] <- list()
}

## use criterion TwiceMean (second choice)
pcdimmat2 <- matrix(0, M, 19) # save pc components
bicdimmat2 <- matrix(0, M, 19) # save number of clusters based on BIC
noiselist2 <- list() # save noise features
for (i in 1:19) {
  noiselist2[[i]] <- list()
}

clustlist2 <- list() # save features' grouping labels
for (i in 1:19) {
  clustlist2[[i]] <- list()
}

## define indices in NbClust function
indfun <- c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", 
           "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", 
           "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", 
           "gamma", "gplus", "tau", "dunn", "sdindex", "sdbw")

# indfun <- c("kl", "ch", "hartigan", "cindex", "db", "silhouette",      
#           "ratkowsky", "ball", "ptbiserial", "gap", "mcclain", 
#           "gamma", "gplus", "tau", "dunn", "sdindex", "sdbw")

# save number of clusters computed from Nbclust indices
res <- matrix(0, 19*M, length(indfun))
cpt.time <- matrix(0, M, 19)
twicemean.time <- matrix(0, M, 19)
scod.time <- matrix(0, M, 19)
nbclust.time <- list()
for (i in 1:length(indfun)) {
  nbclust.time[[i]] <- matrix(0, M, 19)
}


###-----------------------------------------------------------------------------------------
### Simulation (M=1000 datasets for each correlation structure)
### In manuscript, we are using data types 4 - 19

set.seed(123)
# ptm <- proc.time()
for (m in 1:M) {
  cat("Iteration", m, "\n", file=stderr())
  flush.console()
  ## firstly use thresher with CPT (l=7) criterion

  # simulated data type 1
  mu <- rnorm(p, sd = 9)
  sim1 <- mvrnorm(N, mu, cormatlist[[1]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim1)))), agfuns)
  ptm <- proc.time()
  thresher1 <- Thresher(sim1, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper1 <- Reaper(thresher1, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,1] <- as.numeric(t)[3]
  pcdimmat[m,1] <- thresher1@pcdim 
  noiselist[[1]][[m]] <- (1:ncol(sim1))[thresher1@delta<=0.3] 
  bicdimmat[m,1] <- testreaper1@nGroups
  if (!any(is.na(testreaper1@fit))) {
    clustlist[[1]][[m]] <- apply(testreaper1@fit$P, 1, which.max)
  } else {
    clustlist[[1]][[m]] <- NA
  }

  # simulated data type 2
  sim2 <- matrix(rgamma(N*p, shape=2, scale=5), N, p, byrow=FALSE)
  sim2 <- sweep(sim2-10, 2, sqrt(lambdavec), "*")/sqrt(50)
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim2)))), agfuns)
  ptm <- proc.time()
  thresher2 <- Thresher(sim2, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper2 <- Reaper(thresher2, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,2] <- as.numeric(t)[3]
  pcdimmat[m,2] <- thresher2@pcdim 
  noiselist[[2]][[m]] <- (1:ncol(sim2))[thresher2@delta<=0.3]
  bicdimmat[m,2] <- testreaper2@nGroups
  if (!any(is.na(testreaper2@fit))) {
    clustlist[[2]][[m]] <- apply(testreaper2@fit$P, 1, which.max)
  } else {
    clustlist[[2]][[m]] <- NA
  }

  # simulated data type 3
  df <- 3
  sim3 <- matrix(rt(N*p, df), N, p, byrow=FALSE)
  sim3 <- sweep(sim3, 2, sqrt(lambdavec), "*")/sqrt(3)
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim3)))), agfuns)
  ptm <- proc.time()
  thresher3 <- Thresher(sim3, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper3 <- Reaper(thresher3, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,3] <- as.numeric(t)[3]
  pcdimmat[m,3] <- thresher3@pcdim
  noiselist[[3]][[m]] <- (1:ncol(sim3))[thresher3@delta<=0.3] 
  bicdimmat[m,3] <- testreaper3@nGroups
  if (!any(is.na(testreaper3@fit))) {
    clustlist[[3]][[m]] <- apply(testreaper3@fit$P, 1, which.max)
  } else {
    clustlist[[3]][[m]] <- NA
  }

  # simulated data type 4
  mu <- rnorm(p, sd = 9)
  sim4 <- mvrnorm(N, mu, cormatlist[[4]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim4)))), agfuns)
  ptm <- proc.time()
  thresher4 <- Thresher(sim4, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper4 <- Reaper(thresher4, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,4] <- as.numeric(t)[3]
  pcdimmat[m,4] <- thresher4@pcdim
  noiselist[[4]][[m]] <- (1:ncol(sim4))[thresher4@delta<=0.3]
  bicdimmat[m,4] <- testreaper4@nGroups 
  if (!any(is.na(testreaper4@fit))) {
    clustlist[[4]][[m]] <- apply(testreaper4@fit$P, 1, which.max)
  } else {
    clustlist[[4]][[m]] <- NA
  }

  # simulated data type 5
  mu <- rnorm(p, sd = 9)
  sim5 <- mvrnorm(N, mu, cormatlist[[5]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim5)))), agfuns)
  ptm <- proc.time()
  thresher5 <- Thresher(sim5, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper5 <- Reaper(thresher5, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,5] <- as.numeric(t)[3]
  pcdimmat[m,5] <- thresher5@pcdim
  noiselist[[5]][[m]] <- (1:ncol(sim5))[thresher5@delta<=0.3] 
  bicdimmat[m,5] <- testreaper5@nGroups 
  if (!any(is.na(testreaper5@fit))) {
    clustlist[[5]][[m]] <- apply(testreaper5@fit$P, 1, which.max)
  } else {
    clustlist[[5]][[m]] <- NA
  }

  # simulated data type 6
  mu <- rnorm(p, sd = 9)
  sim6 <- mvrnorm(N, mu, cormatlist[[6]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim6)))), agfuns)
  ptm <- proc.time()
  thresher6 <- Thresher(sim6, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper6 <- Reaper(thresher6, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,6] <- as.numeric(t)[3]
  pcdimmat[m,6] <- thresher6@pcdim
  noiselist[[6]][[m]] <- (1:ncol(sim6))[thresher6@delta<=0.3] 
  bicdimmat[m,6] <- testreaper6@nGroups
  if (!any(is.na(testreaper6@fit))) {
    clustlist[[6]][[m]] <- apply(testreaper6@fit$P, 1, which.max)
  } else {
    clustlist[[6]][[m]] <- NA
  }

  # simulated data type 7
  mu <- rnorm(p, sd = 9)
  sim7 <- mvrnorm(N, mu, cormatlist[[7]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim7)))), agfuns)
  ptm <- proc.time()
  thresher7 <- Thresher(sim7, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper7 <- Reaper(thresher7, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,7] <- as.numeric(t)[3]
  pcdimmat[m,7] <- thresher7@pcdim
  noiselist[[7]][[m]] <- (1:ncol(sim7))[thresher7@delta<=0.3]
  bicdimmat[m,7] <- testreaper7@nGroups
  if (!any(is.na(testreaper7@fit))) {
    clustlist[[7]][[m]] <- apply(testreaper7@fit$P, 1, which.max)
  } else {
    clustlist[[7]][[m]] <- NA
  }

  # simulated data type 8
  mu <- rnorm(p, sd = 9)
  sim8 <- mvrnorm(N, mu, cormatlist[[8]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim8)))), agfuns)
  ptm <- proc.time()
  thresher8 <- Thresher(sim8, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper8 <- Reaper(thresher8, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,8] <- as.numeric(t)[3]
  pcdimmat[m,8] <- thresher8@pcdim
  noiselist[[8]][[m]] <- (1:ncol(sim8))[thresher8@delta<=0.3]
  bicdimmat[m,8] <- testreaper8@nGroups 
  if (!any(is.na(testreaper8@fit))) {
    clustlist[[8]][[m]] <- apply(testreaper8@fit$P, 1, which.max)
  } else {
    clustlist[[8]][[m]] <- NA
  }

  # simulated data type 9
  mu <- rnorm(p, sd = 9)
  sim9 <- mvrnorm(N, mu, cormatlist[[9]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim9)))), agfuns)
  ptm <- proc.time()
  thresher9 <- Thresher(sim9, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper9 <- Reaper(thresher9, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,9] <- as.numeric(t)[3]
  pcdimmat[m,9] <- thresher9@pcdim
  noiselist[[9]][[m]] <- (1:ncol(sim9))[thresher9@delta<=0.3] 
  bicdimmat[m,9] <- testreaper9@nGroups
  if (!any(is.na(testreaper9@fit))) {
    clustlist[[9]][[m]] <- apply(testreaper9@fit$P, 1, which.max)
  } else {
    clustlist[[9]][[m]] <- NA
  }

  # simulated data type 10
  mu <- rnorm(p, sd = 9)
  sim10 <- mvrnorm(N, mu, cormatlist[[10]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim10)))), agfuns)
  ptm <- proc.time()
  thresher10 <- Thresher(sim10, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper10 <- Reaper(thresher10, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,10] <- as.numeric(t)[3]
  pcdimmat[m,10] <- thresher10@pcdim
  noiselist[[10]][[m]] <- (1:ncol(sim10))[thresher10@delta<=0.3]
  bicdimmat[m,10] <- testreaper10@nGroups 
  if (!any(is.na(testreaper10@fit))) {
    clustlist[[10]][[m]] <- apply(testreaper10@fit$P, 1, which.max)
  } else {
    clustlist[[10]][[m]] <- NA
  }

  # simulated data type 11
  mu <- rnorm(p, sd = 9)
  sim11 <- mvrnorm(N, mu, cormatlist[[11]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim11)))), agfuns)
  ptm <- proc.time()
  thresher11 <- Thresher(sim11, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper11 <- Reaper(thresher11, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,11] <- as.numeric(t)[3]
  pcdimmat[m,11] <- thresher11@pcdim
  noiselist[[11]][[m]] <- (1:ncol(sim11))[thresher11@delta<=0.3]
  bicdimmat[m,11] <- testreaper11@nGroups 
  if (!any(is.na(testreaper11@fit))) {
    clustlist[[11]][[m]] <- apply(testreaper11@fit$P, 1, which.max)
  } else {
    clustlist[[11]][[m]] <- NA
  }

  # simulated data type 12
  mu <- rnorm(p, sd = 9)
  sim12 <- mvrnorm(N, mu, cormatlist[[12]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim12)))), agfuns)
  ptm <- proc.time()
  thresher12 <- Thresher(sim12, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper12 <- Reaper(thresher12, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,12] <- as.numeric(t)[3]
  pcdimmat[m,12] <- thresher12@pcdim
  noiselist[[12]][[m]] <- (1:ncol(sim12))[thresher12@delta<=0.3]
  bicdimmat[m,12] <- testreaper12@nGroups 
  if (!any(is.na(testreaper12@fit))) {
    clustlist[[12]][[m]] <- apply(testreaper12@fit$P, 1, which.max)
  } else {
    clustlist[[12]][[m]] <- NA
  }

  # simulated data type 13
  mu <- rnorm(p, sd = 9)
  sim13 <- mvrnorm(N, mu, cormatlist[[13]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim13)))), agfuns)
  ptm <- proc.time()
  thresher13 <- Thresher(sim13, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper13 <- Reaper(thresher13, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,13] <- as.numeric(t)[3]
  pcdimmat[m,13] <- thresher13@pcdim
  noiselist[[13]][[m]] <- (1:ncol(sim13))[thresher13@delta<=0.3]
  bicdimmat[m,13] <- testreaper13@nGroups 
  if (!any(is.na(testreaper13@fit))) {
    clustlist[[13]][[m]] <- apply(testreaper13@fit$P, 1, which.max)
  } else {
    clustlist[[13]][[m]] <- NA
  }

  # simulated data type 14
  mu <- rnorm(p, sd = 9)
  sim14 <- mvrnorm(N, mu, cormatlist[[14]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim14)))), agfuns)
  ptm <- proc.time()
  thresher14 <- Thresher(sim14, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper14 <- Reaper(thresher14, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,14] <- as.numeric(t)[3]
  pcdimmat[m,14] <- thresher14@pcdim
  noiselist[[14]][[m]] <- (1:ncol(sim14))[thresher14@delta<=0.3] 
  bicdimmat[m,14] <- testreaper14@nGroups 
  if (!any(is.na(testreaper14@fit))) {
    clustlist[[14]][[m]] <- apply(testreaper14@fit$P, 1, which.max)
  } else {
    clustlist[[14]][[m]] <- NA
  }

  # simulated data type  15
  mu <- rnorm(p, sd = 9)
  sim15 <- mvrnorm(N, mu, cormatlist[[15]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim15)))), agfuns)
  ptm <- proc.time()
  thresher15 <- Thresher(sim15, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper15 <- Reaper(thresher15, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,15] <- as.numeric(t)[3]
  pcdimmat[m,15] <- thresher15@pcdim
  noiselist[[15]][[m]] <- (1:ncol(sim15))[thresher15@delta<=0.3] 
  bicdimmat[m,15] <- testreaper15@nGroups  
  if (!any(is.na(testreaper15@fit))) {
    clustlist[[15]][[m]] <- apply(testreaper15@fit$P, 1, which.max)
  } else {
    clustlist[[15]][[m]] <- NA
  }

  # simulated data type 16
  mu <- rnorm(p, sd = 9)
  sim16 <- mvrnorm(N, mu, cormatlist[[16]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim16)))), agfuns)
  ptm <- proc.time()
  thresher16 <- Thresher(sim16, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper16 <- Reaper(thresher16, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,16] <- as.numeric(t)[3]
  pcdimmat[m,16] <- thresher16@pcdim
  noiselist[[16]][[m]] <- (1:ncol(sim16))[thresher16@delta<=0.3] 
  bicdimmat[m,16] <- testreaper16@nGroups 
  if (!any(is.na(testreaper16@fit))) {
    clustlist[[16]][[m]] <- apply(testreaper16@fit$P, 1, which.max)
  } else {
    clustlist[[16]][[m]] <- NA
  }

  # simulated data type 17
  mu <- rnorm(p, sd = 9)
  sim17 <- mvrnorm(N, mu, cormatlist[[17]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim17)))), agfuns)
  ptm <- proc.time()
  thresher17 <- Thresher(sim17, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper17 <- Reaper(thresher17, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,17] <- as.numeric(t)[3]
  pcdimmat[m,17] <- thresher17@pcdim
  noiselist[[17]][[m]] <- (1:ncol(sim17))[thresher17@delta<=0.3] 
  bicdimmat[m,17] <- testreaper17@nGroups 
  if (!any(is.na(testreaper17@fit))) {
    clustlist[[17]][[m]] <- apply(testreaper17@fit$P, 1, which.max)
  } else {
    clustlist[[17]][[m]] <- NA
  }

  # simulated data type 18
  mu <- rnorm(p, sd = 9)
  sim18 <- mvrnorm(N, mu, cormatlist[[18]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim18)))), agfuns)
  ptm <- proc.time()
  thresher18 <- Thresher(sim18, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper18 <- Reaper(thresher18, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,18] <- as.numeric(t)[3]
  pcdimmat[m,18] <- thresher18@pcdim
  noiselist[[18]][[m]] <- (1:ncol(sim18))[thresher18@delta<=0.3] 
  bicdimmat[m,18] <- testreaper18@nGroups 
  if (!any(is.na(testreaper18@fit))) {
    clustlist[[18]][[m]] <- apply(testreaper18@fit$P, 1, which.max)
  } else {
    clustlist[[18]][[m]] <- NA
  }

  # simulated data type 19
  mu <- rnorm(p, sd = 9)
  sim19 <- mvrnorm(N, mu, cormatlist[[19]])
  #   compareAgDimMethods(AuerGervini(SamplePCA(t(scale(sim19)))), agfuns)
  ptm <- proc.time()
  thresher19 <- Thresher(sim19, method="auer.gervini", scale=TRUE, agfun=agfuns[[l]])
  testreaper19 <- Reaper(thresher19, useLoadings=TRUE)
  t <- proc.time() - ptm
  cpt.time[m,19] <- as.numeric(t)[3]
  pcdimmat[m,19] <- thresher19@pcdim
  noiselist[[19]][[m]] <- (1:ncol(sim19))[thresher19@delta<=0.3] 
  bicdimmat[m,19] <- testreaper19@nGroups 
  if (!any(is.na(testreaper19@fit))) {
    clustlist[[19]][[m]] <- apply(testreaper19@fit$P, 1, which.max)
  } else {
    clustlist[[19]][[m]] <- NA
  }


  ## secondly use thresher with TwiceMean criterion

  ptm <- proc.time()
  thresher_v1 <- Thresher(sim1, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v1 <- Reaper(thresher_v1, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,1] <- as.numeric(t)[3]
  pcdimmat2[m,1] <- thresher_v1@pcdim 
  noiselist2[[1]][[m]] <- (1:ncol(sim1))[thresher_v1@delta<=0.3] 
  bicdimmat2[m,1] <- testreaper_v1@nGroups
  if (!any(is.na(testreaper_v1@fit))) {
    clustlist2[[1]][[m]] <- apply(testreaper_v1@fit$P, 1, which.max)
  } else {
    clustlist2[[1]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v2 <- Thresher(sim2, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v2 <- Reaper(thresher_v2, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,2] <- as.numeric(t)[3]
  pcdimmat2[m,2] <- thresher_v2@pcdim 
  noiselist2[[2]][[m]] <- (1:ncol(sim2))[thresher_v2@delta<=0.3]
  bicdimmat2[m,2] <- testreaper_v2@nGroups
  if (!any(is.na(testreaper_v2@fit))) {
    clustlist2[[2]][[m]] <- apply(testreaper_v2@fit$P, 1, which.max)
  } else {
    clustlist2[[2]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v3 <- Thresher(sim3, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v3 <- Reaper(thresher_v3, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,3] <- as.numeric(t)[3]
  pcdimmat2[m,3] <- thresher_v3@pcdim
  noiselist2[[3]][[m]] <- (1:ncol(sim3))[thresher_v3@delta<=0.3] 
  bicdimmat2[m,3] <- testreaper_v3@nGroups
  if (!any(is.na(testreaper_v3@fit))) {
    clustlist2[[3]][[m]] <- apply(testreaper_v3@fit$P, 1, which.max)
  } else {
    clustlist2[[3]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v4 <- Thresher(sim4, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v4 <- Reaper(thresher_v4, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,4] <- as.numeric(t)[3]
  pcdimmat2[m,4] <- thresher_v4@pcdim
  noiselist2[[4]][[m]] <- (1:ncol(sim4))[thresher_v4@delta<=0.3]
  bicdimmat2[m,4] <- testreaper_v4@nGroups 
  if (!any(is.na(testreaper_v4@fit))) {
    clustlist2[[4]][[m]] <- apply(testreaper_v4@fit$P, 1, which.max)
  } else {
    clustlist2[[4]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v5 <- Thresher(sim5, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v5 <- Reaper(thresher_v5, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,5] <- as.numeric(t)[3]
  pcdimmat2[m,5] <- thresher_v5@pcdim
  noiselist2[[5]][[m]] <- (1:ncol(sim5))[thresher_v5@delta<=0.3] 
  bicdimmat2[m,5] <- testreaper_v5@nGroups 
  if (!any(is.na(testreaper_v5@fit))) {
    clustlist2[[5]][[m]] <- apply(testreaper_v5@fit$P, 1, which.max)
  } else {
    clustlist2[[5]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v6 <- Thresher(sim6, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v6 <- Reaper(thresher_v6, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,6] <- as.numeric(t)[3]
  pcdimmat2[m,6] <- thresher_v6@pcdim
  noiselist2[[6]][[m]] <- (1:ncol(sim6))[thresher_v6@delta<=0.3] 
  bicdimmat2[m,6] <- testreaper_v6@nGroups
  if (!any(is.na(testreaper_v6@fit))) {
    clustlist2[[6]][[m]] <- apply(testreaper_v6@fit$P, 1, which.max)
  } else {
    clustlist2[[6]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v7 <- Thresher(sim7, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v7 <- Reaper(thresher_v7, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,7] <- as.numeric(t)[3]
  pcdimmat2[m,7] <- thresher_v7@pcdim
  noiselist2[[7]][[m]] <- (1:ncol(sim7))[thresher_v7@delta<=0.3]
  bicdimmat2[m,7] <- testreaper_v7@nGroups 
  if (!any(is.na(testreaper_v7@fit))) {
    clustlist2[[7]][[m]] <- apply(testreaper_v7@fit$P, 1, which.max)
  } else {
    clustlist2[[7]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v8 <- Thresher(sim8, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v8 <- Reaper(thresher_v8, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,8] <- as.numeric(t)[3]
  pcdimmat2[m,8] <- thresher_v8@pcdim
  noiselist2[[8]][[m]] <- (1:ncol(sim8))[thresher_v8@delta<=0.3]
  bicdimmat2[m,8] <- testreaper_v8@nGroups 
  if (!any(is.na(testreaper_v8@fit))) {
    clustlist2[[8]][[m]] <- apply(testreaper_v8@fit$P, 1, which.max)
  } else {
    clustlist2[[8]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v9 <- Thresher(sim9, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v9 <- Reaper(thresher_v9, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,9] <- as.numeric(t)[3]
  pcdimmat2[m,9] <- thresher_v9@pcdim
  noiselist2[[9]][[m]] <- (1:ncol(sim9))[thresher_v9@delta<=0.3] 
  bicdimmat2[m,9] <- testreaper_v9@nGroups
  if (!any(is.na(testreaper_v9@fit))) {
    clustlist2[[9]][[m]] <- apply(testreaper_v9@fit$P, 1, which.max)
  } else {
    clustlist2[[9]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v10 <- Thresher(sim10, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v10 <- Reaper(thresher_v10, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,10] <- as.numeric(t)[3]
  pcdimmat2[m,10] <- thresher_v10@pcdim
  noiselist2[[10]][[m]] <- (1:ncol(sim10))[thresher_v10@delta<=0.3]
  bicdimmat2[m,10] <- testreaper_v10@nGroups 
  if (!any(is.na(testreaper_v10@fit))) {
    clustlist2[[10]][[m]] <- apply(testreaper_v10@fit$P, 1, which.max)
  } else {
    clustlist2[[10]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v11 <- Thresher(sim11, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v11 <- Reaper(thresher_v11, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,11] <- as.numeric(t)[3]
  pcdimmat2[m,11] <- thresher_v11@pcdim
  noiselist2[[11]][[m]] <- (1:ncol(sim11))[thresher_v11@delta<=0.3]
  bicdimmat2[m,11] <- testreaper_v11@nGroups 
  if (!any(is.na(testreaper_v11@fit))) {
    clustlist2[[11]][[m]] <- apply(testreaper_v11@fit$P, 1, which.max)
  } else {
    clustlist2[[11]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v12 <- Thresher(sim12, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v12 <- Reaper(thresher_v12, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,12] <- as.numeric(t)[3]
  pcdimmat2[m,12] <- thresher_v12@pcdim
  noiselist2[[12]][[m]] <- (1:ncol(sim12))[thresher_v12@delta<=0.3]
  bicdimmat2[m,12] <- testreaper_v12@nGroups
  if (!any(is.na(testreaper_v12@fit))) {
    clustlist2[[12]][[m]] <- apply(testreaper_v12@fit$P, 1, which.max)
  } else {
    clustlist2[[12]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v13 <- Thresher(sim13, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v13 <- Reaper(thresher_v13, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,13] <- as.numeric(t)[3]
  pcdimmat2[m,13] <- thresher_v13@pcdim
  noiselist2[[13]][[m]] <- (1:ncol(sim13))[thresher_v13@delta<=0.3]
  bicdimmat2[m,13] <- testreaper_v13@nGroups 
  if (!any(is.na(testreaper_v13@fit))) {
    clustlist2[[13]][[m]] <- apply(testreaper_v13@fit$P, 1, which.max)
  } else {
    clustlist2[[13]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v14 <- Thresher(sim14, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v14 <- Reaper(thresher_v14, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,14] <- as.numeric(t)[3]
  pcdimmat2[m,14] <- thresher_v14@pcdim
  noiselist2[[14]][[m]] <- (1:ncol(sim14))[thresher_v14@delta<=0.3] 
  bicdimmat2[m,14] <- testreaper_v14@nGroups 
  if (!any(is.na(testreaper_v14@fit))) {
    clustlist2[[14]][[m]] <- apply(testreaper_v14@fit$P, 1, which.max)
  } else {
    clustlist2[[14]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v15 <- Thresher(sim15, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v15 <- Reaper(thresher_v15, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,15] <- as.numeric(t)[3]
  pcdimmat2[m,15] <- thresher_v15@pcdim
  noiselist2[[15]][[m]] <- (1:ncol(sim15))[thresher_v15@delta<=0.3] 
  bicdimmat2[m,15] <- testreaper_v15@nGroups 
  if (!any(is.na(testreaper_v15@fit))) {
    clustlist2[[15]][[m]] <- apply(testreaper_v15@fit$P, 1, which.max)
  } else {
    clustlist2[[15]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v16 <- Thresher(sim16, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v16 <- Reaper(thresher_v16, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,16] <- as.numeric(t)[3]
  pcdimmat2[m,16] <- thresher_v16@pcdim
  noiselist2[[16]][[m]] <- (1:ncol(sim16))[thresher_v16@delta<=0.3] 
  bicdimmat2[m,16] <- testreaper_v16@nGroups
  if (!any(is.na(testreaper_v16@fit))) {
    clustlist2[[16]][[m]] <- apply(testreaper_v16@fit$P, 1, which.max)
  } else {
    clustlist2[[16]][[m]] <- NA
  }  

  ptm <- proc.time()
  thresher_v17 <- Thresher(sim17, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v17 <- Reaper(thresher_v17, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,17] <- as.numeric(t)[3]
  pcdimmat2[m,17] <- thresher_v17@pcdim
  noiselist2[[17]][[m]] <- (1:ncol(sim17))[thresher_v17@delta<=0.3] 
  bicdimmat2[m,17] <- testreaper_v17@nGroups 
  if (!any(is.na(testreaper_v17@fit))) {
    clustlist2[[17]][[m]] <- apply(testreaper_v17@fit$P, 1, which.max)
  } else {
    clustlist2[[17]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v18 <- Thresher(sim18, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v18 <- Reaper(thresher_v18, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,18] <- as.numeric(t)[3]
  pcdimmat2[m,18] <- thresher_v18@pcdim
  noiselist2[[18]][[m]] <- (1:ncol(sim18))[thresher_v18@delta<=0.3] 
  bicdimmat2[m,18] <- testreaper_v18@nGroups  
  if (!any(is.na(testreaper_v18@fit))) {
    clustlist2[[18]][[m]] <- apply(testreaper_v18@fit$P, 1, which.max)
  } else {
    clustlist2[[18]][[m]] <- NA
  }

  ptm <- proc.time()
  thresher_v19 <- Thresher(sim19, method="auer.gervini", scale=TRUE, agfun=agfuns[[1]])
  testreaper_v19 <- Reaper(thresher_v19, useLoadings=TRUE)
  t <- proc.time() - ptm
  twicemean.time[m,19] <- as.numeric(t)[3]
  pcdimmat2[m,19] <- thresher_v19@pcdim
  noiselist2[[19]][[m]] <- (1:ncol(sim19))[thresher_v19@delta<=0.3] 
  bicdimmat2[m,19] <- testreaper_v19@nGroups 
  if (!any(is.na(testreaper_v19@fit))) {
    clustlist2[[19]][[m]] <- apply(testreaper_v19@fit$P, 1, which.max)
  } else {
    clustlist2[[19]][[m]] <- NA
  }


  ## compute number of clusters based on NbClust indices

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim1), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,1] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+1,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+1,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim2), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,2] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+2,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+2,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim3), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,3] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+3,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+3,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim4), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,4] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+4,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+4,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim5), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,5] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+5,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+5,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim6), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,6] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+6,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+6,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim7), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,7] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+7,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+7,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim8), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,8] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+8,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+8,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim9), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,9] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+9,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+9,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim10), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,10] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+10,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+10,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim11), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,11] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+11,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+11,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim12), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,12] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+12,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+12,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim13), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,13] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+13,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+13,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim14), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,14] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+14,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+14,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim15), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,15] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+15,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+15,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim16), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,16] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+16,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+16,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim17), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,17] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+17,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+17,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim18), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,18] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+18,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+18,i] <- NA }
  }

  for (i in 1:length(indfun)) {
    ptm <- proc.time()
    fit <- try(NbClust(t(sim19), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
            method = "ward.D2", index = indfun[i]))
    t <- proc.time() - ptm
    nbclust.time[[i]][m,19] <- as.numeric(t)[3]
    if (class(fit)!="try-error") {
    res[19*(m-1)+19,i] <- fit$Best.nc[1] 
    } else { res[19*(m-1)+19,i] <- NA }
  }

  ## finally use SCOD algorithm to get clusters and outliers

  ptm <- proc.time()
  scod_1 <- SCOD(sim1)
  t <- proc.time() - ptm
  scod.noiselist[[1]][[m]] <- scod_1$outliers 
  scod.ncmat[m,1] <- scod_1$nclust
  scod.clustlist[[1]][[m]] <- scod_1$clasI * (scod_1$clasI <= scod_1$nclust)
  scod.time[m,1] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_2 <- SCOD(sim2)
  t <- proc.time() - ptm
  scod.noiselist[[2]][[m]] <- scod_2$outliers 
  scod.ncmat[m,2] <- scod_2$nclust
  scod.clustlist[[2]][[m]] <- scod_2$clasI * (scod_2$clasI <= scod_2$nclust)
  scod.time[m,2] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_3 <- SCOD(sim3)
  t <- proc.time() - ptm
  scod.noiselist[[3]][[m]] <- scod_3$outliers 
  scod.ncmat[m,3] <- scod_3$nclust
  scod.clustlist[[3]][[m]] <- scod_3$clasI * (scod_3$clasI <= scod_3$nclust)
  scod.time[m,3] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_4 <- SCOD(sim4)
  t <- proc.time() - ptm
  scod.noiselist[[4]][[m]] <- scod_4$outliers 
  scod.ncmat[m,4] <- scod_4$nclust
  scod.clustlist[[4]][[m]] <- scod_4$clasI * (scod_4$clasI <= scod_4$nclust)
  scod.time[m,4] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_5 <- SCOD(sim5)
  t <- proc.time() - ptm
  scod.noiselist[[5]][[m]] <- scod_5$outliers 
  scod.ncmat[m,5] <- scod_5$nclust
  scod.clustlist[[5]][[m]] <- scod_5$clasI * (scod_5$clasI <= scod_5$nclust)
  scod.time[m,5] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_6 <- SCOD(sim6)
  t <- proc.time() - ptm
  scod.noiselist[[6]][[m]] <- scod_6$outliers 
  scod.ncmat[m,6] <- scod_6$nclust
  scod.clustlist[[6]][[m]] <- scod_6$clasI * (scod_6$clasI <= scod_6$nclust)
  scod.time[m,6] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_7 <- SCOD(sim7)
  t <- proc.time() - ptm
  scod.noiselist[[7]][[m]] <- scod_7$outliers 
  scod.ncmat[m,7] <- scod_7$nclust
  scod.clustlist[[7]][[m]] <- scod_7$clasI * (scod_7$clasI <= scod_7$nclust)
  scod.time[m,7] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_8 <- SCOD(sim8)
  t <- proc.time() - ptm
  scod.noiselist[[8]][[m]] <- scod_8$outliers 
  scod.ncmat[m,8] <- scod_8$nclust
  scod.clustlist[[8]][[m]] <- scod_8$clasI * (scod_8$clasI <= scod_8$nclust)
  scod.time[m,8] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_9 <- SCOD(sim9)
  t <- proc.time() - ptm
  scod.noiselist[[9]][[m]] <- scod_9$outliers 
  scod.ncmat[m,9] <- scod_9$nclust
  scod.clustlist[[9]][[m]] <- scod_9$clasI * (scod_9$clasI <= scod_9$nclust)
  scod.time[m,9] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_10 <- SCOD(sim10)
  t <- proc.time() - ptm
  scod.noiselist[[10]][[m]] <- scod_10$outliers 
  scod.ncmat[m,10] <- scod_10$nclust
  scod.clustlist[[10]][[m]] <- scod_10$clasI * (scod_10$clasI <= scod_10$nclust)
  scod.time[m,10] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_11 <- SCOD(sim11)
  t <- proc.time() - ptm
  scod.noiselist[[11]][[m]] <- scod_11$outliers 
  scod.ncmat[m,11] <- scod_11$nclust
  scod.clustlist[[11]][[m]] <- scod_11$clasI * (scod_11$clasI <= scod_11$nclust)
  scod.time[m,11] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_12 <- SCOD(sim12)
  t <- proc.time() - ptm
  scod.noiselist[[12]][[m]] <- scod_12$outliers 
  scod.ncmat[m,12] <- scod_12$nclust
  scod.clustlist[[12]][[m]] <- scod_12$clasI * (scod_12$clasI <= scod_12$nclust)
  scod.time[m,12] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_13 <- SCOD(sim13)
  t <- proc.time() - ptm
  scod.noiselist[[13]][[m]] <- scod_13$outliers 
  scod.ncmat[m,13] <- scod_13$nclust
  scod.clustlist[[13]][[m]] <- scod_13$clasI * (scod_13$clasI <= scod_13$nclust)
  scod.time[m,13] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_14 <- SCOD(sim14)
  t <- proc.time() - ptm
  scod.noiselist[[14]][[m]] <- scod_14$outliers 
  scod.ncmat[m,14] <- scod_14$nclust
  scod.clustlist[[14]][[m]] <- scod_14$clasI * (scod_14$clasI <= scod_14$nclust)
  scod.time[m,14] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_15 <- SCOD(sim15)
  t <- proc.time() - ptm
  scod.noiselist[[15]][[m]] <- scod_15$outliers 
  scod.ncmat[m,15] <- scod_15$nclust
  scod.clustlist[[15]][[m]] <- scod_15$clasI * (scod_15$clasI <= scod_15$nclust)
  scod.time[m,15] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_16 <- SCOD(sim16)
  t <- proc.time() - ptm
  scod.noiselist[[16]][[m]] <- scod_16$outliers 
  scod.ncmat[m,16] <- scod_16$nclust
  scod.clustlist[[16]][[m]] <- scod_16$clasI * (scod_16$clasI <= scod_16$nclust)
  scod.time[m,16] <- as.numeric(t)[3] 

  ptm <- proc.time()
  scod_17 <- SCOD(sim17)
  t <- proc.time() - ptm
  scod.noiselist[[17]][[m]] <- scod_17$outliers 
  scod.ncmat[m,17] <- scod_17$nclust
  scod.clustlist[[17]][[m]] <- scod_17$clasI * (scod_17$clasI <= scod_17$nclust)
  scod.time[m,17] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_18 <- SCOD(sim18)
  t <- proc.time() - ptm
  scod.noiselist[[18]][[m]] <- scod_18$outliers 
  scod.ncmat[m,18] <- scod_18$nclust
  scod.clustlist[[18]][[m]] <- scod_18$clasI * (scod_18$clasI <= scod_18$nclust)
  scod.time[m,18] <- as.numeric(t)[3]

  ptm <- proc.time()
  scod_19 <- SCOD(sim19)
  t <- proc.time() - ptm
  scod.noiselist[[19]][[m]] <- scod_19$outliers 
  scod.ncmat[m,19] <- scod_19$nclust
  scod.clustlist[[19]][[m]] <- scod_19$clasI * (scod_19$clasI <= scod_19$nclust)
  scod.time[m,19] <- as.numeric(t)[3]

}
# proc.time() - ptm

# estimated computation time = 11 hours

# assign column names for results matrix
colnames(res) <- indfun

# save results in workspace 
save.image("simulation1.RData")




