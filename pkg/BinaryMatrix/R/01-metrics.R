### Definitions of various similarity or distance measures between
### binary vectors. In the main reeference (below), Choi and colleagues
### reviewed 76 different measures of simialrity or distance between
### binary vectors. They also produced a hierarchical clustering of
### these measures, based on the correlation between their distance
### values on multiple simulated data sets. For metrics that are highly
### similar, we choose a single representative.


# [1] Choi SS, Cha SH, Tappert CC, A Survey of Binary Similarity and Distance Measures.
#     Systemics, Cybernetics, and Informatics. 2010; 8(1):43-48.

### Cluster 1 contains Dice & Sorenson, Ochiai, Kulczynski, Bray & Curtis,
### Baroni-Urbani & Buser, and Jaccard 
jaccardSimilarity <- function(X) {
  N11 <- X %*% t(X)
  N01 <- (1 - X) %*% t(X)
  N10 <- X %*% t(1 - X)
  den <- (N11 + N01 + N10)
  den[den == 0] <- 1
  N11 / den
}
jaccardDistance <- function(X) {
  as.dist( 1 - jaccardSimilarity(X))
}

### Cluster 2 contains Sokal & Sneath, Gilbert & Wells, Gower & Legendre,
### Pearson & Heron, and Sokal & Michener.
sokalMichenerSimilarity <- function(X) {
  N11 <- X %*% t(X)
  N01 <- (1 - X) %*% t(X)
  N10 <- X %*% t(1 - X)
  N00 <- (1 - X) %*% t(1 - X)
  (N00 + N11) / (N00 + N11 + N01 + N10)
}
sokalMichenerDistance <- function(X) {3
  as.dist( 1 - sokalMichenerSimilarity(X))
}

### Also in Cluster 2 are Hamming, Manhattan, Canberra, Minkowski,
### and Euclidean. We might want to try Euclidean as well.
### Hamming
hammingDistance <- function(X) {
  b <- (1 - X) %*% t(X)
  c <- X %*% t(1 - X)
  as.dist(b + c)
}

### Cluster 3 contains Driver & Kroeber, Forbes, Fossum, and
### Russell & Rao
russellRaoSimilarity <- function(X) {
  N11 <- X %*% t(X)
  N11/ ncol(X)
}
russellRaoDistance <- function(X) {
  as.dist( 1 - russellRaoSimilarity(X) )
}

### Remaining metrics are more isolated, without strong clutering.
### We consider a few examples.

### Pearson correlation
pearsonSimilarity <- function(X) {
  a <- X %*% t(X)
  b <- (1 - X) %*% t(X)
  c <- X %*% t(1 - X)
  d <- (1 - X) %*% t(1 - X)
  n <- a + b + c + d
  (n*(a*d - b*c)^2)/((a + b)*(a + c)*(b + d)*(c + d))
}
pearsonDistance <- function(X) {
  as.dist( ncol(X) - pearsonSimilarity(X) )
}

### Goodman & Kruskal
goodmanKruskalSimilarity <- function(X) {
  a <- X %*% t(X)
  b <- (1 - X) %*% t(X)
  c <- X %*% t(1 - X)
  d <- (1 - X) %*% t(1 - X)
  n <- a + b + c + d
  sig <- pmax(a, b) + pmax(c, d) + pmax(a, c) + pmax(b, d)
  sigprime <- pmax(a + c) + pmax(a + b, c + d)
  (sig - sigprime) / (2*n - sigprime)
}
goodmanKruskalDistance <- function(X) {
  as.dist( 1 - goodmanKruskalSimilarity(X) )
}

binaryDistance <- function(X, metric) {
  METRICS <- c(jaccard = jaccardDistance,
               sokalMichener = sokalMichenerDistance,
               hamming = hammingDistance,
               russellRao = russellRaoDistance,
               pearson = pearsonDistance,
               goodmanKruskal = goodmanKruskalDistance,
               manhattan = function(x) {
                    DD <- dist(t(X), "manhattan")
                    DD/max(DD)
               },
               canberra = function(x) dist(t(X), "canberra")/ncol(X),
               binary = function(x) dist(t(X), "binary"),
               euclid = function(x) {
                    DD <- dist(t(X), "euclid")
                    DD/max(DD)
               }
               ) # why is this a list and not a vector?
  fullname <- match.arg(metric, names(METRICS))
  FUN <- METRICS[[fullname]]
  RES <- FUN(t(X))
  comment(RES) <- fullname
  RES
}
