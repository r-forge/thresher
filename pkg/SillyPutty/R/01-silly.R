setOldClass("silhouette")

setClass("SillyPutty",
         slots = c(
           cluster = "numeric",
           silhouette = "silhouette",
           minSW = "numeric",
           meanSW = "numeric")
)

SillyPutty <- function(labels, dissim, maxIter = 1000, loopSize = 15, verbose = FALSE) {
  refined <- labels
  sily <- silhouette(refined, dist = dissim)
  temp <- as.matrix(sily)
  counter <- 0
  minis <- NULL
  asw <- NULL
  repeat {
    x <- temp[w <- which.min(temp[,3]),]
    minis <- c(minis, x[3])
    refinednew <- refined
    refinednew[w] <- x[2]
    silynew <- silhouette(refinednew, dist = dissim)
    asw <- mean(silynew[,3])
    tempnew <- as.matrix(silynew)
    if(verbose) cat("Prev:\t", x[3], "\tCurr:\t", min(tempnew[,3]), "\n", file=stderr())
    refined <- refinednew
    temp <- tempnew
    counter <-  counter + 1
    if (min(tempnew[,3]) >= 0) break    # all well classified
    L <- length(minis)
    if (L > loopSize) {
      if (min(tempnew[,3]) %in% minis[(L - loopSize):L]) {
        minis <- c(minis, min(tempnew[,3]))
        gack <<- list(refined=refined, minis=minis)
        warning("Early halt because of possible infinite loop.")
        break
      }
    }
    if (counter > maxIter) {
      warning("maximum number of iterations reached.")
      break;
    }
  }
  new("SillyPutty",
      cluster = refined,
      silhouette = silynew,
      minSW = minis,
      meanSW = asw)
}

setClass("RandomSillyPutty",
         slots = c(
           MX = "numeric",
           MN = "numeric",
           ave = "numeric",
           labels = "matrix",
           minItemSW = "list"
         ))

# add a "verbose" argument, so you can tell if there is any progress
RandomSillyPutty <- function(distobj, K, N = 100, verbose = FALSE, ...) {
  res <- list()
  for (tries in 1:N) {
    if (verbose) cat(tries, ' ', file=stderr())
    START <- sample(K, attr(distobj, "Size"), replace=TRUE)
    FINAL <- SillyPutty(START, distobj, ...)
    if (length(unique(FINAL@cluster)) < K) next
    res[[1 + length(res)]] <- FINAL
  }
  if (verbose) cat("\n", file=stderr())
  ave <- sapply(res, function(R) R@meanSW)
  mx <- which.max(ave)[1]
  mn <- which.min(ave)[1]
  MN <- res[[mn]]@cluster
  MX <- res[[mx]]@cluster
  MN <- remap(MX, MN)
  labs <- lapply(res, function(R) R@cluster)
  labs <- t(as.matrix(as.data.frame(lapply(labs,
                                          function(L) remap(MX, L)))))
  colnames(labs) <- attr(distobj, "Labels")
  new("RandomSillyPutty",
      MX = MX, MN = MN, ave = ave,
      labels = labs,
      minItemSW = lapply(res, function(R) R@minSW))
}

summarizeRepeats <- function(multiple) {
  lab <- multiple@labels
  ave <- multiple@ave[which(!duplicated(lab))]
  u <- unique(lab)
  rownames(u) <- 1:nrow(u)

  weight <- apply(u, 1, function(v) sum(apply(lab, 1, function(x) all(x == v))))
  common <- u[which.max(weight)[1],]
  rsc <- cyanyellow(16)[cut(ave, 16, labels=FALSE)]
  rsc[which.max(weight)] <- "red"
  list(MODE = common, MIN = multiple@MN, MAX=multiple@MX,
       u = u, ave = ave, lab = lab, weight = weight, rsc = rsc)
}

showme <- function(clue1, clue2, visual, distobj, palette = NULL) {
  NG1 <- length(unique(clue1))
  NG2 <- length(unique(clue2))
  M <- max(NG1, NG2)
  if (is.null(palette)) {
    if (M > 44) {
      stop("Unable to generate a palette with more than 44 colors.")
    } else if (M > 36) {
      D <- dark.colors(24)
      L <- light.colors(24)
      palette <- sortByHue(c(D,L))
      names(palette) <- colorNames(palette)
      palette <-palette[!duplicated(palette)]
    } else {
      palette <- palette36.colors(M)
    }
  } else {
    if (M > length(palette)) {
      stop("Palette does not contain enough colors.")
    } # else, use what you've got
  }
  n1 <- deparse(substitute(clue1))
  n2 <- deparse(substitute(clue2))
  opar <- par(mfcol=c(2, 2))
  on.exit(par(opar))
  plot(visual, pch=16, col=palette[clue1], cex=1.7, main=n1)
  plot(silhouette(clue1, dist = distobj), col=palette[1:NG1])
  plot(visual, pch=16, col=palette[clue2], cex=1.7, main=n2)
  plot(silhouette(clue2, dist = distobj), col=palette[1:NG2])
}

