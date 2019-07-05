### Based on the code from SPADE by Peng Qiu, which was implemented
### specifically for masa cyteometry data stored in FCS files.
downsample <- function(target, distanceMat, cutoff) {
  localDensity <- function(distanceMat, epsilon) {
    apply(distanceMat, 1, function(x) mean(x < epsilon))
  }
  ldens <- localDensity(distanceMat, cutoff)
  scalefactor <- sum(1/ldens)/length(ldens)
  prob <- (1/ldens)/length(ldens) * target/scalefactor
  prob[prob > 0.99999] <- 0.99999
  P <- rbinom(length(ldens), 1, prob)
  P == 1
}

createGraph <- function(DV, Q) {
  M <- as.matrix(DV@distance) # spread out distance matrix
  U <- M[upper.tri(M)] # then reselect the upper triangular portion
  if (missing(Q)) Q <- quantile(M, 0.1)
  selU <- U < Q        # select edges based on a distance cutoff
  ## Compute column indices for selected edges
  M <- 0*M
  M <- sweep(M, 2, 1:ncol(M), "+")
  cols <- M[upper.tri(M)][selU]
  ## doi the same for row indices
  M <- 0*M
  M <- sweep(M, 1, 1:ncol(M), "+")
  rows <- M[upper.tri(M)][selU]
  ## then build an edge-list data frame
  daft <- data.frame(A = colnames(M)[cols], 
                     B = colnames(M)[rows], 
                     weight = 1-U[selU])
  ## make an igraph from the edges
  myg <- graph_from_data_frame(daft, directed=FALSE)
  myg <- set_vertex_attr(myg, "size", value=3)   # shrink the nodes
  myg <- set_vertex_attr(myg, "label", value="") # hide the labels
  V <- vertex_attr(myg)
#  syms <- c("square", "circle")[DV@symv - 14] # BUG!
  syms <- c("square", "circle")[1 + (DV@symv %% 2)]
  names(syms) <- names(DV@symv)
  myg <- set_vertex_attr(myg, "color", value=DV@colv[V$name])
  myg <- set_vertex_attr(myg, "shape", value=syms[V$name])
  layouts <- list(nicely = layout_nicely(myg))
  if (!is.null(MV <- DV@view[["mds"]])) {
    V <- igraph::vertex_attr(myg)
    layouts[["mds"]] <- MV[V$name,]
  }
  if (!is.null(TV <- DV@view[["tsne"]])) {
    Y <- TV$Y
    rownames(Y) <- labels(DV@distance)
    V <- igraph::vertex_attr(myg)
    layouts[["tsne"]] <- Y[V$name,]
  }
  list(graph = myg, layouts = layouts)
}
