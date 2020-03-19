setClass("Mercator",
         slots = c(
           metric ="character",
           distance = "dist",
           view = "list",
           palette = "character",
           symbols = "numeric",
           clusters = "numeric")
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
                        # since it removes one leaf at a tim
    as.hclust(P)
  }
  switch(name,
         mds = object[i,],
         tsne = shrinkTSNE(object, i),
         hclust = snipTree(object, i),
         heat = object
         )
}

setMethod("[", signature = "Mercator", function(x, i, j, ..., drop=FALSE) {
  # ignore j and drop, since we are thinking of this as really one-dimenaional
  if(missing(i)) i <- j
  M <- as.matrix(x@distance)[i,i]
  V <- list()
  for (I in 1:length(x@view)) {
    N <- names(x@view)[I]
    cat(N, "\n", file=stderr())
    V[[N]] <- shrinkView(N, x@view[[I]], i)
  }
  new("Mercator",
      metric = x@metric,
      distance = as.dist(M),
      view = V,
      palette = x@palette,
      symbols = x@symbols,
      clusters = x@clusters[i]
      )
})

setMethod("dim", signature = "Mercator", function(x) {
  dim(as.matrix((x@distance)))
})

setMethod("summary", signature(object="Mercator"),
          function(object, ...) {
  cat("An object of the 'Mercator' class, using the '", object@metric, "' metric, of size\n")
  print(dim(object))
  cat("Contains these visualizations: ", names(object@view), "\n")
})

setMethod("hist", signature(x = "Mercator"), function(x, breaks=123, ...) {
  M <- as.matrix(x@distance)
  U <- M[upper.tri(M)]
  hist(U, breaks=breaks, ...)
})

setMethod("barplot", signature("Mercator"),
          function(height, main = '', sub = NULL, border = NA, space = 0, ...) {
  clus <- getClusters(height)
  silh <- silhouette(clus, height@distance)
  od <- order(clus, silh[,'sil_width'])
  if (is.null(sub)) {
    swd <- silh[od, 'sil_width']
    sub <- paste("Mean SW =", round(mean(swd), 5))
  }
  barplot(swd, col = colv(height)[od], border=NA, 
          sub=sub, main=main, ...)
})

setMethod("plot", signature("Mercator", "missing"),
          function(x, view = NULL, ask = NULL, ...) {
### known kinds of visualizations
  plotMDS <- function(x, ...) {
    plot(x@view[["mds"]], col=colv(x), pch=symv(x), xlab = "PC1", ylab="PC2", ...)
  }
  plotTSNE <- function(x, ...) {
    plot(x@view[["tsne"]]$Y, col=colv(x), pch=symv(x), xlab = "T1", ylab="T2", ...)
  }
  plotHC <- function(x, ...) {
    dend <- x@view[["hclust"]]
    labs <- dend$labels
    if (is.null(labs)) {
      labs <- seq(1, length(dend$order))
    }
    plotColoredClusters(dend, cols = colv(x), labs = labs, ...)
  }
  plotIG <- function(x, layout = NULL, ...) {
    G <- x@view[["graph"]]
    if (is.null(layout)) {
      layout <- G$layouts[[1]]
    }
    if (is.character(layout)) {
      foo <- ifelse(layout %in% c("mds", "tsne"),
                    jitter,
                    function(x) x)
      layout <- foo(G$layouts[[layout[1]]])
    }
    plot(G$graph, layout = layout, ...)
  }
### implications of 'view' and 'ask' parameters
  if (is.null(view)) { # first attached view is the default
    view <- list(names(x@view)[1])
  }
  if (length(view) ==1 & view == "all") {
    view <- names(x@view)
  }
  if (!is.list(view)) { # can show more than one
    view <- as.list(view)
  }
  if (is.null(ask)) {
    ask <- prod(par("mfcol")) < length(view) && dev.interactive()
  }
  if (ask & length(view) > 1) { # ask politely
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
### actually plot stuff
  for (V in view) {
    switch(V,
           mds = plotMDS(x, ...),
           tsne = plotTSNE(x, ...),
           hclust = plotHC(x, ...),
           graph = plotIG(x, ...))
  }
  invisible(x)
})

setMethod("scatter", signature(object = "Mercator"),
          function(object, view = NULL, ask = NULL,
                   colramp = NULL, ...) {
### known kinds of visualizations
  smoothMDS <- function(object, ...) {
    smoothScatter(x = object@view[["mds"]], xlab = "PC1", ylab="PC2", ...)
  }
  smoothTSNE <- function(object, ...) {
    smoothScatter(x = object@view[["tsne"]]$Y, xlab = "T1", ylab="T2", ...)
  }
### implications of 'view' and 'ask' parameters
  if (is.null(view)) { # first attached view is the default
    view <- list(names(object@view)[1])
  }
  if (length(view) ==1 & view == "all") {
    view <- names(object@view)
  }
  if (!is.list(view)) { # can show more than one
    view <- as.list(view)
  }
  if (is.null(ask)) {
    ask <- prod(par("mfcol")) < length(view) && dev.interactive()
  }
  if (ask && length(view) > 1) { # ask politely
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
### actually smooth stuff
  if (is.null(colramp)) {
    colramp <- topo.colors
  }
  for (V in view) {
    switch(V,
           mds = smoothMDS(object, colramp = colramp, ...),
           tsne = smoothTSNE(object, colramp = colramp, ...),
           cat("No smooth scatter plot is available for view '", V, "'.\n"))
  }
  invisible(object)
})


makeDisplay <- function(object) {
  pal  <-  object@palette
  baseSyms  <-  object@symbols
  clusters <- object@clusters
  if (is.null(names(clusters))) {
    names(clusters) <- attr(object@distance, "Labels")
  }
  K <- max(clusters)
  L <- length(pal)
  R <- ifelse(K %% L == 0, K/L, 1 + trunc(K/L))
  if (R > length(baseSyms)) {
    stop("Are you kidding me? You can't possibly want that many (", K, ") clusters.")
  }
  mycol <- rep(pal, R)
  mysym <- rep(baseSyms[1:R], each=L)
  colv <- mycol[clusters]
  symv <- mysym[clusters]
  names(colv) <- names(symv) <- names(clusters)
  list(colv = colv, symv = symv)
}

RC <- function(colv, symv) {
  Dark24 <- dark.colors(24)
  baseSyms <- c(16, 15, 17, 18, 10, 7, 11, 9)
  lead <- sapply(symv, function(X) which(baseSyms == X)) - 1
  units <- sapply(colv, function(X) which(Dark24 == X))
  24*lead + units
}

colv <- function(object) {
  md <- makeDisplay(object)
  md$colv
}
symv <- function(object) {
  md <- makeDisplay(object)
  md$symv
}

getClusters <- function(DV) {
  DV@clusters
}

remapColors <- function(fix, vary) {
  fixCluster <- getClusters(fix)
  varyCluster <- getClusters(vary)
  newCluster <- remap(fixCluster, varyCluster)
  names(newCluster) <- names(varyCluster)
  new("Mercator",
      metric = vary@metric,
      distance = vary@distance,
      view = vary@view,
      palette = vary@palette,
      symbols = vary@symbols,
      clusters = newCluster
      )
}

setClusters <- function(DV, clusters) {
  if (is.null(names(clusters))) {
    names(clusters) <- names(DV@clusters)
  }
  L <- attr(DV@distance, "Labels")
  if (any(names(clusters) != L)) {
    if (all(names(clusters) %in% L)) {
      warning("Reordering cluster info to match existing names.")
      clusters <- clusters[L]
    } else {
      stop("Cluster names do not match existing names.")
    }
  }
  new("Mercator",
      metric = DV@metric,
      distance = DV@distance,
      view = DV@view,
      palette = DV@palette,
      symbols = DV@symbols,
      clusters = clusters
      )
}
recolor <- function(DV, clusters) {
  .Deprecated("setClusters",
              msg = "Please use the new, more descriptive name (setClusters) for this function.")
  setClusters(DV, clusters)
}

Mercator <- function(binaryMat, metric, method, K, ...) {
  DistMat <- binaryDistance(binaryMat@binmat, metric)
  M <- attr(DistMat, "comment")
  if (!is.null(M)) metric <- M
  clust <- pam(DistMat, k=K, diss=TRUE, cluster.only=TRUE)
  if (is.null(names(clust))) names(clust) <- attr(DistMat, "Labels")
  view <- list()
  ob <- new("Mercator",
            metric = metric,
            distance = DistMat,
            view = view,
            palette = dark.colors(24),
            symbols = c(16, 15, 17, 18, 10, 7, 11, 9),
            clusters = clust
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
  temp <- do.call(FUN, argList)
  if (is.matrix(temp) && is.null(rownames(temp))) {
    rownames(temp) <- attr(DV@distance, "Labels")
  }
  DV@view[[method]] <- temp
  DV
}

