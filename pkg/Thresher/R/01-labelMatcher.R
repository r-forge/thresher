labelMatcher <- function(tab, verbose=FALSE) {
  if (ncol(tab) != nrow(tab)) stop("must supply a square matrix or table")
  N <- nrow(tab)
  if (N==1) return(list(ii=1, jj=1))
  if (all(tab==0)) return(list(ii=1:N, jj=1:N))

  rowm <- apply(tab, 1, function(x) {
    js <- which(x==max(x))
    sap <- sapply(js, function(j, x) {
      y <- tab[,j]
      tot <- sum(x) + sum(y) - x[j]
      x[j]/tot
    }, x=x)
    max(sap)
  })
  i0 <- which(rowm==max(rowm, na.rm=TRUE))[1]
  j0 <- which(tab[i0,]==max(tab[i0,]))[1]
  mat <- tab[-i0, -j0, drop=FALSE]
  if (verbose) {
    cat("i0 =", i0, "j0 =", j0, "\n")
    print(dim(mat))
    print(mat)
  }
  recurse <- labelMatcher(mat)
  idx <- (1:N)[-i0][recurse$ii]
  jdx <- (1:N)[-j0][recurse$jj]
  value <- list(ii=c(i0, idx), jj=c(j0, jdx))
  value
}

matchLabels <- function(tab) {
  matches <- labelMatcher(tab)
  tab[matches$ii, matches$jj]
}

countAgreement <- function(tab) {
  sum(diag( matchLabels(tab) ))
}

labelAccuracy <- function(data, labels, linkage="ward.D2") {
  # order matters: we pick the first of 'most accurate' as best
  metrics <- c("pearson", "spearman", "euclidean",
               "uncentered correlation",  "absolute pearson",
               "sqrt pearson",
               "minkowski", "maximum",  "weird")
  labels <- as.factor(labels)
  nGroups <- length(levels(labels))
  accu <- sapply(metrics, function(m) {
    hc <- hclust(distanceMatrix(data, m, p=1), linkage)
    tab <- table(labels, paste("K", cutree(hc, k=nGroups), sep=""))  
    countAgreement(tab)
  })
  accu/ncol(data)
}

bestMetric <- function(data, labels) {
  accu <- labelAccuracy(data, labels)
  names(accu)[accu==max(accu)][1]
}

remap <- function(fix, vary) {
  if(is.factor(vary)) { # need consistent names
    V <- vary
  } else {
    V <- factor(vary)
  }
  tab <- table(fix, vary)
  lem <- labelMatcher(tab)
  oj <- order(lem$jj)
  tricky <- lem$ii[oj]
  names(tricky) <- levels(V)
  res <- levels(V)[tricky[vary]] # should swap the labels in "vary" to best match "fix"
  if (is.numeric(vary)) res <- as.numeric(res)
  res
}
