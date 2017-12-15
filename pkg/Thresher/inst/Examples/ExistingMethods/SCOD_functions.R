
### functions to implement the SCOD algorithm

## ModifySimilarity - function to modify the similarity considering the top 15 % 
## similarities S is the original similarity

ModifySimilarity <- function(S) {
  
  N <- nrow(S)
  K <- ceiling(0.15 * N) 
  A <- apply(S, 2, sort)
  W <- colSums(A[(N-K+1):N, ])/K
  # MinMax normalize the weights
  WT <- W
  W <- (WT - min(WT))/(max(WT) - min(WT))
  WW <- rep(1, N) %*% t(W) 
  S1 <- S * WW * t(WW)
  S1 <- (S1 - min(S1))/(max(S1) - min(S1))

  return(S1)

}

## replicator - Replicator dynamics to return the equilibrium point (solution to opt.)
## Parameters: A - N X N similarity matrix; x - N X 1 starting point in the simplex;
## toll - stop when the velocity of the dynamics is below toll;
## maxiter - (optional) stop after maxiter iterations; 
## Returned values: x - N X 1 point in the simplex representing the cluster; count - 
## number of required iterations.

replicator <- function(A, x, toll = 1e-6, maxiter = 1000) {

  error <- 2 * toll + 1

  count <- 0
  while(error > toll & count < maxiter) {
      x_old <- x
      x <- x * (A %*% x)
      if (sum(x) > 0) {
        x <- x/sum(x)
      } else {
        x <- x_old
        break
      }
      error <- sqrt(sum((x - x_old)^2))
      count <- count + 1
  }
  
  return(list(x = x, count = count))

}


## SCOD - SCOD algorithm to get the number of clusters and outliers
## data is the dataset that need to be clustered, columns are objects
SCOD <- function(data, sigma = 0.8, toll = 1e-5, maxiter = 1000) {

  D <- as.matrix(dist(t(data), diag = TRUE))
  D <- D/max(D)
  A <- exp(-(D/sigma)^2)
  A <- (A - min(A))/(max(A) - min(A))
  diag(A) <- 0
  
  n <- nrow(A)
  x <- runif(n) + 10000000
  x <- x/sum(x)
  S <- ModifySimilarity(A)
  GC <- t(x) %*% A %*% x
  idx <- 1:n
  clasI <- idx
  iter <- 1
  cohesiVandSize <- NULL
  
  nclust <- 0
  outlier.set <- NULL
  cluster.set <- list()

  while (length(idx) > 1) {
    
    x <- runif(n) + 10000000
    x <- x/sum(x)
    replicator.res <- replicator(A, x, toll, maxiter)
    y <- replicator.res[[1]]
    C <- t(y) %*% A %*% y
    CL <- t(y) %*% S %*% y

    cluster <- which(y > 1e-5)

    if (C < GC | CL < GC) {
      outlier.set <- c(outlier.set, idx[cluster])      
    } else {
      nclust <- nclust + 1
      cluster.set[[nclust]] <- idx[cluster]
    }

    MC <- A[cluster, cluster]
    diag(MC) <- max(MC)
    MinS <- min(MC)
    Lnew <- length(cluster)
    SL <- nrow(A)
    if (length(cluster) < SL - 1) {
      A <- A[-cluster, ]
      A <- A[ , -cluster]    
      S <- S[-cluster, ]
      S <- S[ , -cluster]
    } else { 
      A <- matrix(A[-cluster, ], nrow = SL - length(cluster), ncol = SL)
      A <- matrix(A[ , -cluster], nrow = SL - length(cluster), ncol = SL - length(cluster))
      S <- matrix(S[-cluster, ], nrow = SL - length(cluster), ncol = SL)
      S <- matrix(S[ , -cluster], nrow = SL - length(cluster), ncol = SL - length(cluster))
    }    

    # MaxS <- max(A)
    n <- length(idx[-cluster])
    cohesiVandSize <- rbind(cohesiVandSize, t(c(iter, C, CL)))
    clasI[idx[cluster]] <- iter
    idx <- idx[-cluster]
    iter <- iter + 1

  }

  if (length(idx) == 1) {
      outlier.set <- c(outlier.set, idx)
  }

  return(list(clasI = clasI, cohesiVandSize = cohesiVandSize, nclust = nclust, clusters = 
         cluster.set, outliers = outlier.set))

}



