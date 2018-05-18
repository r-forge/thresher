setOldClass("movMF")
setClassUnion("fit or miss", c("movMF", "logical"))
setClassUnion("number or miss", c("numeric", "logical"))

.fitModels <- function(thresher, useLoadings) {
  D <- thresher@pcdim
  top <- max(2*D+1, 2)
  bot <- ifelse(D==0, 1, D)
  ncrange <- bot:top 
  fits <- list()
  for (nc in ncrange) {
    K <- ifelse(D > 1, D, 2)
    if (useLoadings) {
      top <- ifelse(nc > 2*K, K+1, K)
      space <- thresher@loadings[, 1:top]
    } else {
      space <- t(scale(thresher@data))
    }
    fitted <- try(movMF(space, nc, start="S"), silent=TRUE)
    if (inherits(fitted, "try-error")) {
      fitted <- try(movMF(space, nc, runs=10), silent=TRUE)
      if (inherits(fitted, "try-error")) {
        next
      }
    }
    M <- 1+length(fits)
    fits[[M]] <- fitted
    names(fits)[M] <- paste("NC=", nc, sep="")
  }
  fits
}


setClass("SignalSet",
         slots = c(members="list",
                   continuous="matrix",
                   binary="matrix",
                   continuousClusters="hclust",
                   binaryClusters="hclust"
                   ))

.findSignals <- function(thresher, fit, nGroups, linkage="ward.D2") {
  gassign <- predict(fit) # groups of features defined by vMF mixture.
  # get feature names in each group
  featureNames <- lapply(1:nGroups, function(idx) {
    colnames(thresher@data)[gassign==idx]
  })
  # compute the average across each group in the space of loadings
  foo <- as.matrix(as.data.frame(lapply(1:nGroups, function(idx) {
    y <- thresher@loadings[gassign==idx, 1:thresher@pcdim, drop=FALSE]
    apply(y, 2, mean)
  })))
  if (any(is.nan(foo))) {
    cat("Robot FU!!\n", file = stderr())
    for (ii in 1:nrow(foo)) {
      cat(paste(foo[ii, ], collapse = "\t"), "\n", file = stderr())
    }
    cat("gassign:", paste(gassign, collapse = " "), "\n", 
        file = stderr())
    cat("pcdim:", thresher@pcdim, "\n", file = stderr())
    stop("Danger, Will Robinson")
  }
  # order signals by strength along PCs
  if (ncol(foo) > 1) {
    temp <- sqrt(apply(foo^2, 2, cumsum))
    if (class(temp) == "numeric") temp <- matrix(temp, nrow=1)
    cumulative <- apply(temp, 2, sum)
    oc <- rev(order(cumulative)) # biggest is first
    featureNames <- featureNames[oc]
    foo <- foo[, oc, drop=FALSE]
  }
  colnames(foo) <- paste("G", 1:nGroups, sep="")
  # convert to directions
  u <- unitize(foo)
  mat <- t(u)%*%u
  # account for round-off error that can blowup acos
  mat[mat < -1] <- -1
  mat[mat > 1] <- 1
  degrees <- round(acos(mat)*360/(2*pi))
  diag(degrees) <- 0

  # figure out which signals point in opposite directions
  pcSamples <- thresher@spca@scores[, 1:thresher@pcdim]
  mark <- rep(FALSE, nGroups)
  signals <- list()
  members <- list()
  while(TRUE) {
    if(all(mark)) break
    w <- which(!mark)[1]
    oppose <- (!mark & degrees[w,] > 135)
    if (any(oppose)) {
      v <- which(oppose)
      if (length(v) > 1) {
        warning("should not be opposite to two different things")
        md <- max(degrees[w,][!mark])
        v <- which(!mark & degrees[w,]==md)[1]
      }
      signal <- matrix(foo[,w] - foo[,v], ncol=1)
      meme <- list(LEFT=featureNames[w], RIGHT=featureNames[v])
      mark[c(v, w)] <- TRUE
    } else {
      signal <- matrix(foo[,w], ncol=1)
      meme <- list(FEAT=featureNames[w])
      mark[w] <- TRUE
    }
    M <- length(signals)
    signals[[1+M]] <- pcSamples %*% signal
    members[[1+M]] <- meme
  }
  contSignal <- as.matrix(as.data.frame(signals))
  colnames(contSignal) <- paste("Sig", 1:length(signals), sep="")
  conClust <- hclust(distanceMatrix(t(contSignal), "euclid"), linkage)
  binSignal <- 1*(contSignal > 0)
  binClust <- hclust(dist.binary(binSignal, method=2), linkage)
  new("SignalSet",
      members=members,
      continuous=contSignal,
      binary=binSignal,
      continuousClusters=conClust,
      binaryClusters=binClust)
}

setClass("Reaper",
         contains = "Thresher",
         slots = c(useLoadings="logical",
                   keep="logical",
                   nGroups="number or miss",
                   fit="fit or miss",
                   allfits="list",
                   bic="number or miss",
                   metric="character",
                   signalSet="SignalSet",
                   maxSampleGroups="numeric"
                   ))

Reaper <- function(thresher,
                   useLoadings=FALSE,
                   cutoff=0.3,
                   metric=NULL,
                   linkage="ward.D2",
                   maxSampleGroups=0,
                   ...) {
  keep <- thresher@delta > cutoff
  m <- ifelse(is.null(metric), "pearson", metric)
  cleaned <- Thresher(thresher@data[, keep],
                      paste(thresher@name, "cleaned", sep='.'),
                      metric=m,
                      method="auer.gervini",
                      agfun=thresher@agfun,
                      ...)
  tab <- 0
  counter <- 0
  while(any(tab==0) & counter < 5) {
    counter <- counter + 1
    fits <- .fitModels(cleaned, useLoadings)
    if (length(fits)==0)  next
    bic <- sapply(fits, BIC)
    woo <- which(bic==min(bic))
    ng <- as.integer(sub("NC=", "", names(woo)))
    fit <- fits[[woo]]
    # sometimes there are no samples assigned to one of the groups,
    # which causes later steps to fail catastrophically
    gassign <- factor(predict(fit), levels=1:ng)
    tab <- table(gassign)
  }
  if (length(fits)==0) {
    bic <- ng <- fit <- NA
    metric <- "no fit"
    sigset <- new("SignalSet")
  } else {
    if (is.null(metric)) {
      pp <- factor(paste("G", predict(fit), sep=""))
      metric <- bestMetric(cleaned@data, pp)
      cleaned@gc <- hclust(distanceMatrix(cleaned@data, metric, p=1), linkage)
    }
    if (any(tab==0)) {
      sigset <- new("SignalSet")
    } else {
      sigset <- .findSignals(cleaned, fit, ng, linkage=linkage)
    }
  }
  new("Reaper", cleaned, useLoadings=useLoadings,
      keep=keep, nGroups=ng, fit=fit,
      allfits=fits, bic=bic, metric=metric,
      signalSet=sigset, maxSampleGroups=maxSampleGroups)
}

# adds to REVERSE thresherPalette to define SAMPLE colors
.makeColorScheme <- function(N) {
  grays <- paste("gray", 1:100, sep="")
  colorscheme <- c(samplePalette, grays)[1:N]
  if (any(is.na(colorscheme))) {
    n <- sum(is.na(colorscheme))
    fresh <- colors()
    fresh <- fresh[!(fresh %in%  grays)]
    colorscheme[is.na(colorscheme)] <- sample(fresh, n)
  }
  colorscheme
}

# defines SAMPLE colors
setMethod("getSplit", "Reaper", function(object) {
  binSignal <- object@signalSet@binary
  contSignal <- object@signalSet@continuous
  nSig <- ncol(binSignal)
  powers <- 2^(-1 + nSig:1)
  if (object@maxSampleGroups > 0) {
    top <- trunc(log2(object@maxSampleGroups))
    if (nSig >= top+1) {
      L <- length(powers)
      powers <- 2^((top-1):0)
      powers <- c(powers, rep(0, L - length(powers)))
    }
  }
  weights <- matrix(powers, ncol=1)
  sclass <- as.vector(binSignal %*% weights)
  op <- order(sclass, contSignal[,1])
  sclass <- factor(paste("C", sclass))
  colorscheme <- .makeColorScheme(length(levels(sclass)))
  colset <- colorscheme[as.numeric(sclass)]
  fc <- factor(colset, levels=colorscheme)
  list(fc=fc, op=op, colset=colset)
})

# really? we either select groups from hierarchical clustering OR
# from the vMF fits, without making it explicit?

# uses thresherPalette to define PROTEIN/ANTIBODY colors
setMethod("getColors", "Reaper", function(object, K=NULL) {
  if(is.null(K)) {
    K <- object@nGroups
    abnum <- predict(object@fit)
  } else {
    abnum <- cutree(object@gc, k=K) # do we need this here?
  }
  L <- length(thresherPalette)
  if (K > L) {
    abnum <- abnum %% L
    abnum[abnum==0] <- L
  }
  thresherPalette[abnum]
})

setMethod("getStyles", "Reaper", function(object, K=NULL) {
  if(is.null(K)) K <- object@nGroups
  abnum <- predict(object@fit)
  1 +  trunc(abnum/length(thresherPalette))
})

unitize <- function(mat) {
  enorm <- sqrt(apply(mat^2, 2, sum))
  sweep(mat, 2, enorm, "/")
}

setMethod("makeFigures", "Reaper", function(object, DIR=NULL, main=NULL, ...) {
  if(is.null(main)) main <- object@name
  callNextMethod()
  if (is.null(DIR)) opar <- par(ask=TRUE)
  fname <- gsub("\\.", "-", object@name) # latex-safe
  gs <- getSplit(object)
  K <- length(levels(gs$fc))
  colsch <- .makeColorScheme(K)
# fig6
  if (!is.null(DIR)) {
    png(filename=file.path(DIR, paste(fname, "-06-heat-cont.png", sep="")),
        width=1930, height=640, pointsize=12, bg="white")
  }
  hc <- object@signalSet@continuousClusters
  csc <- colsch[cutree(hc, k=K)]
  heat(object, Colv=as.dendrogram(hc), ColSideColors=csc,
       main=paste(main, "(Continuous Signals)"),
       wExp=3, margins=c(5,15),
       ...)
  if(!is.null(DIR)) dev.off()
# fig7
  if (!is.null(DIR)) {
    png(filename=file.path(DIR, paste(fname, "-07-heat-bin.png", sep="")),
        width=1920, height=640, pointsize=12, bg="white")
  }
  hc <- object@signalSet@binaryClusters
  csc <- colsch[cutree(hc, k=K)]
  heat(object, Colv=as.dendrogram(hc), ColSideColors=csc,
       main=paste(main, "(Binary Signals)"),
       wExp=3, margins=c(5,15),
       ...)
  if(!is.null(DIR)) dev.off()
  if (is.null(DIR)) par(opar)
  invisible(object)
})
