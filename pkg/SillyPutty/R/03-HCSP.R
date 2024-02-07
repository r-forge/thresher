### Copyright (C) Kevin R. Coombes, 2024

HCSP <- function(dis, K, method = "ward.D2", ...) {
  if (!inherits(dis, "dist")) {
    dis <- dist(dis) # defaults to Euclidean
  }
  hc <- hclust(dis, method = method)
  sp <- SillyPutty(cutree(hc, k = K, ...), dis)
  list(hc = hc, sp = sp)
}
