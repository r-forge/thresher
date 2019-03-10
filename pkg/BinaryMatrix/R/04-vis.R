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

shrinkView <- function(name, object, i) {
  shrinkTSNE <- function(tis, i) {
    tis$Y <- tis$Y[i,]
    tis
  }
  snipTree <- function(hc, i) {
    dend <- as.dendrogram(hc)
    L <- labels(dend)
    if (is.character(i)) {
      i <- L %in% i
    } else if (is.numeric(i)) {
      temp <- i
      i <- rep(FALSE, length(i))
      i[temp] <- TRUE
    }
    L <- L[!i]
    P <- prune(dend, L) # this step does not scale well
                        # since it removes one leaf at a time
    as.hclust(P)
  }
  switch(name,
         mds = object[i,],
         tsne = shrinkTSNE(object, i),
         hclust = snipTree(object, i),
         heat = object
         )
}

setMethod("[", signature = "DistanceVis", function(x, i, j, ..., drop=FALSE) {
  # ignore j and drop, since we are thinking of this as really one-dimenaional
  if(missing(i)) i <- j
  M <- as.matrix(x@distance)[i,i]
  V <- list()
  for (I in 1:length(x@view)) {
    N <- names(x@view)[I]
    cat(N, "\n", file=stderr())
    V[[N]] <- shrinkView(N, x@view[[I]], i)
  }
  new("DistanceVis",
      metric = x@metric,
      distance = as.dist(M),
      view = V,
      colv = x@colv[i],
      symv = x@symv[i]
      )
})

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
  Dark24 <- dark.colors(24)
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
  Dark24 <- dark.colors(24)
  baseSyms <- c(16, 15, 17, 18, 10, 7, 11, 9)
  lead <- sapply(symv, function(X) which(baseSyms == X)) - 1
  units <- sapply(colv, function(X) which(Dark24 == X))
  24*lead + units
}

remapColors <- function(fix, vary) {
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

recolor <- function(DV, clusters) {
  dispSet <- makeDisplay(clusters)
  colv <- dispSet$colv
  symv <- dispSet$symv
  names(colv) <- names(symv) <- colnames(DV@distance)
  new("DistanceVis",
      metric = DV@metric,
      distance = DV@distance,
      view = DV@view,
      colv = colv,
      symv = symv
      )
}

DistanceVis <- function(binaryMat, metric, method, K, ...) {
  DistMat <- binaryDistance(binaryMat@binmat, metric)
  M <- attr(DistMat, "comment")
  if (!is.null(M)) metric <- M
  poof <- pam(DistMat, k=K, diss=TRUE, cluster.only=TRUE)
  dispSet <- makeDisplay(poof)
  colv <- dispSet$colv
  symv <- dispSet$symv
  names(colv) <- names(symv) <- labels(DistMat)
  view <- list()
  ob <- new("DistanceVis",
            metric = metric,
            distance = DistMat,
            view = view,
            colv = colv,
            symv = symv
            )
  addVisualization(ob, method, ...)
}

addVisualization <- function(DV, method, ...) {
  METHODS = c(mds = function(X, ...) cmdscale(X@distance, ...),
              hclust = function(X, ...) hclust(X@distance, method="ward.D2"),
              tsne = function(X, ...) Rtsne(X@distance, is_distance = TRUE, ...),
              graph = function(X, ...) createGraph(X, ...),
              heat = function(X, ...) X@distance)
  method  <- match.arg(method, names(METHODS))
  if (!is.null(DV@view[[method]])) {
    warning("Overwriting an existing visualization:", method)
  }
  FUN <- METHODS[[method]]
  argList <- c(list(DV), list(...))
  DV@view[[method]] <- do.call(FUN, argList)
  DV
}

