findClusterNumber <- function(distobj, start, end, N = 100,
                              method = c("SillyPutty", "HCSP"), ...) {
  method = match.arg(method)
  res <- list()
  if (method == "SillyPutty") {
    for (k in start:end){
      resMX <- RandomSillyPutty(distobj, k, N, ...)
      tmp <- max(resMX@ave)
      res <- append(res,list(tmp))
    }
    names(res) <- start:end
  } else { # method == "HCSP"
    hc <- hclust(distobj)
    for (k in start:end){
      resMX <- SillyPutty(cutree(hc, k=k), distobj, ...)
      S <- summary(resMX@silhouette)
      tmp <- S$avg.width
      res <- append(res, list(tmp))
    }
    names(res) <- start:end
  }
  res
}
