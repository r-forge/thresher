## D = distance matrix
createSOM <- function(D, grid = somgrid(8, 6, "hexagonal"), ...) {
  Q <- makeEuclidean(D)
  som(Q, grid, ...)
}

createUMAP <- function(D, ...) {
  Q <- makeEuclidean(D)
  umap(Q, ...)
}

makeEuclidean <- function(D) {
  M <- as.matrix(D)
  X <- scale(t(scale(t(M^2), scale = FALSE)), scale = FALSE)
  E <- eigen(-X/2, symmetric = TRUE)$values
  R <- sum(E > 1e-10)
  cmdscale(D, k = R)
}
