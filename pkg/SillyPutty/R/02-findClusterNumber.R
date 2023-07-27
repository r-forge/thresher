findClusterNumber <- function(distobj, start, end, N = 100) {
  set.seed(110)
  res <- list()
  for (k in start:end){
    resMX <- RandomSillyPutty(distobj,k,N)
    tmp <- max(resMX@ave)
    res <- append(res,list(tmp))
  }
  names(res) <- start:end
  res
}
