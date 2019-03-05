library(Rtsne)
library(igraph)
library(Polychrome)

setClass("DistanceVis",
         slots = c(
           metric ="character",
           distance = "dist",
           view = "list",
           colv = "character",
           symv = "numeric")
         )

setMethod("dim", signature = "DistanceVis", function(x) {
  dim(as.matrix((x@distance)))
})

setMethod("summary", signature(object="DistanceVis"),
          function(object, ...) {
  cat("An object of the 'DistanceVis' class, using the '", object@metric, "' metric, of size\n")
  print(dim(object))
  cat("Contains these visualizations: ", names(object@view), "\n")
})

setMethod("hist", signature(x = "DistanceVis"), function(x, breaks=123, ...) {
  M <- as.matrix(x@distance)
  U <- M[upper.tri(M)]
  hist(U, breaks=breaks, ...)
})

makeDisplay <- function(clusters, master =  NULL) {
  data(Dark24, package="Polychrome")
  K <- max(clusters)
  if (!is.null(master)) {
    L <- max(master)
    if (K != L) {
      warning("Mismatch in number of clusters; ignoring.")
    } else {
      clusters <- remap(master, clusters)
    }
  }
  R <- ifelse(K %% 24 == 0, K/24, 1 + trunc(K/24))
  baseSyms <- c(16, 15, 17, 18, 10, 7, 11, 9)
  if (R > length(baseSyms)) {
    stop("Are you kidding me? You can't possibly want that many (", K, ") clusters.")
  }
  mycol <- rep(Dark24, R)
  mysym <- rep(baseSyms[1:R], each=24)
  colv <- mycol[clusters]
  symv <- mysym[clusters]
  list(colv = colv, symv = symv)
}

recoverCluster <- function(colv, symv) {
  baseSyms <- c(16, 15, 17, 18, 10, 7, 11, 9)
  lead <- sapply(symv, function(X) which(baseSyms == X)) - 1
  units <- sapply(colv, function(X) which(Dark24 == X))
  24*lead + units
}

recolor <- function(fix, vary) {
  fixCluster <- recoverCluster(fix@colv, fix@symv)
  varyCluster <- recoverCluster(vary@colv, vary@symv)
  newCluster <- remap(fixCluster, varyCluster)
  newDisplay <- makeDisplay(newCluster)
  new("DistanceVis",
      metric = vary@metric,
      distance = vary@distance,
      view = vary@view,
      colv = newDisplay$colv,
      symv = newDisplay$symv
      )
}

DistanceVis <- function(binaryMat, metric, method, K, ...) {
  DistMat <- binaryDistance(binaryMat@binmat, metric)
  M <- attr(DistMat, "comment")
  if (!is.null(M)) metric <- M
  view <- list()
  view[[method]] <- computeVisualization(DistMat, method, ...)
  poof <- pam(DistMat, k=K, diss=TRUE, cluster.only=TRUE)
  dispSet <- makeDisplay(poof)
  colv <- dispSet$colv
  symv <- dispSet$symv
  names(colv) <- names(symv) <- colnames(DistMat)
  new("DistanceVis",
      metric = metric,
      distance = DistMat,
      view = view,
      colv = colv,
      symv = symv
      )
}

addVisualization <- function(BV, method, ...) {
  if (!is.null(BV@view[[method]])) {
    warning("Overwriting an exisintg visualization:", method)
  }
  BV@view[[method]] <- computeVisualization(BV@distance, method, ...)
  BV
}

computeVisualization <- function(distMat, method, ...) {
  METHODS = c(mds = function(X, ...) cmdscale(X, ...),
              hclust = function(X, ...) hclust(X, method="ward.D2"),
              tsne = function(X, ...) Rtsne(X, is_distance = TRUE, ...),
              heat = function(X, ...) X)
  method  <- match.arg(method, names(METHODS))
  FUN <- METHODS[[method]]
  argList <- c(list(distMat), list(...))
  do.call(FUN, argList)
}
