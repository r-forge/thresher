######################################################
# INTERNAL

spand <- function(v) {
    r <- range(v)
    span <- 0.05*diff(r)
    r1 <- min(0, r[1]-span)
    r2 <- max(r[2]+span, 0)
    c(r1, r2)
}

######################################################
# EXTERNAL

library(colorspace)
N <- 9
foo <- rainbow_hcl(N-1, c=100, l=65)
foo2 <- rainbow(N, s=1, v=1)
baz <- c(foo, foo2, "#aaaaaa", "#eeee00")
thresherPalette <- baz[c(9, 15, 4, 16, 10,
                         18, 19, 13, 8, 6,
                         11, 3, 7, 2, 1,
                         14, 12, 17, 5) ]
rm(N, foo, foo2, baz)
samplePalette <- c("gray", "blue", "red", "purple",
          "green", "cyan", "yellow", "black",
          "gold", "magenta", "brown", "darkcyan",
          "darkgoldenrod", "forestgreen",  "gainsboro", 
          "ghostwhite", "grey100", "honeydew", "ivory4",
          "khaki")

######################################################
# S4 interface

setClass("Thresher",
         slots = c(name="character",
                   data="matrix",
                   spca="SamplePCA",
                   loadings="matrix",
                   gc="hclust",
                   pcdim="numeric",
                   delta="numeric",
                   ag="AuerGervini",
                   agfun='function'
                   ))

Thresher <- function(data,
                     nm=deparse(substitute(data)),
                     metric="pearson",
                     linkage="ward.D2",
                     method=c("auer.gervini", "broken.stick"),
                     scale=TRUE,
                     agfun=agDimTwiceMean) {
  # always center; by default, also 'scale' to standardize
  std <- scale(data, scale=scale)
  spca <- SamplePCA(t(std))
  ag <- AuerGervini(spca)
  method <- match.arg(method)
  pcdim <- switch(method,
                  broken.stick=bsDimension(spca),
                  auer.gervini=agDimension(ag, agfun))
  deltaDim <- max(1, pcdim)
  # distance from origin using loadings
  lambda <- sqrt(spca@variances)
  loadings <- sweep(spca@components, 2, lambda, "*")
  delta  <- sqrt(apply(loadings[,1:deltaDim,drop=FALSE]^2, 1, sum))
  # clustering so we only do it once -- more or less
  gc <- hclust(distanceMatrix(std, metric), linkage)
  new("Thresher",
      name=nm, data=data, spca=spca, loadings=loadings,
      gc=gc, pcdim=pcdim, delta=delta, ag=ag, agfun=agfun)
}

setMethod("screeplot", "Thresher", function(x, col="skyblue", lcol="blue", ...) {
  N <- ncol(x@data)
  bs  <- brokenStick(1:N, N)
  pts <- screeplot(x@spca, main=x@name,
                   xlab="Principal Component",
                   ylab="Fraction of Variance Explained", col=col, ...)
  M <- min(N, length(pts))
  bs <- bs[1:M]
  lines(pts, bs, type='b', lwd=3, pch=16, col=lcol)
  invisible(pts)
})

if (!isGeneric("getColors"))
  setGeneric("getColors",
             function(object, ...) { standardGeneric("getColors") }
             )

# uses thresherPalette to define PROTEIN/ANTIBODY colors
setMethod("getColors", "Thresher", function(object, K=NULL) {
  if (is.null(K)) K <- 1+trunc(log2(ncol(object@data)))
  thresherPalette[cutree(object@gc, k=K)]
})

if (!isGeneric("getStyles"))
  setGeneric("getStyles",
             function(object, ...) { standardGeneric("getStyles") }
             )

setMethod("getStyles", "Thresher", function(object, K=NULL) {
  if (is.null(K)) K <- 1 + trunc(log2(ncol(object@data)))
  gps <- cutree(object@gc, k=K) 
  L <- length(thresherPalette)
  1 + trunc(gps/L)
})

setMethod("plot", c("Thresher", "missing"), function(x, y, ij=1:2, ...) {
  if (length(ij) < 2) ij <- 1:2 # TODO: something better?
  abcols <- getColors(x)
  abstys <- getStyles(x)
  xx <- x@loadings[,ij[1]]
  yy <- x@loadings[,ij[2]]
  plot(xx, yy, type='n', main=x@name,
       xlim=spand(xx), ylim=spand(yy),
       xlab=paste("PC", ij[1], sep=""),
         ylab=paste("PC", ij[2], sp=""))
  for (kdx in 1:nrow(x@loadings)) {
    lines(c(0, xx[kdx]), c(0, yy[kdx]), lwd=3,
          col=abcols[kdx], lty=abstys[kdx]) # here be monsters
  }
  text(xx, yy, rownames(x@loadings), col='black')
  invisible(x)
})

if (!isGeneric("getSplit"))
  setGeneric("getSplit",
             function(object, ...) { standardGeneric("getSplit") }
             )

# use REVERSE thresher palette to define SAMPLE colors
setMethod("getSplit", "Thresher", function(object) {
  colors <- samplePalette[1:4]
  std <- scale(object@data)
  bb <- cutree(object@gc, k=2)
  b1 <- apply(std[,bb==1, drop=FALSE], 1, mean)
  b2 <- apply(std[,bb==2, drop=FALSE], 1, mean)
  c1 <- 1*(b1 > 0)
  c2 <- 1*(b2 > 0)
  colset <- colors[1 + c1 + 2*c2]
  fc <- factor(colset, levels=colors)
  op <- order(c1 + 2*c2)
  list(fc=fc, op=op, colset=colset)
})

if (!isGeneric("scatter"))
  setGeneric("scatter",
             function(object, ...) { standardGeneric("scatter") }
             )

setMethod("scatter", "Thresher", function(object, ...) {
  fc <- getSplit(object)$fc
  plot(object@spca, split=fc, col=levels(fc), main=object@name, ...)
  invisible(object)
})

if (!isGeneric("heat"))
  setGeneric("heat",
             function(object, ...) { standardGeneric("heat") }
             )

# might want to set 'margins=c(5,15)' in ... for long feature names
setMethod("heat", "Thresher", function(object, ncol=NULL,
                                       Colv=NA, main=object@name,
                                       ColSideColors=NULL,
                                       RowSideColors=NULL,
                                       ...) {
  stuff <- getSplit(object)
  if(is.null(ColSideColors)) {
    ColSideColors <- stuff$colset[stuff$op]
  }
  if(is.null(RowSideColors)) {
    RowSideColors <- getColors(object, K=ncol)
  }
  std <- scale(object@data)
  clipped <- std
  K <- 3
  clipped[clipped > K] <- K
  clipped[clipped < -K] <- -K
  forshow <- t(clipped)
  if (any(is.na(Colv))) {
    forshow <- t(clipped[stuff$op,])
  }

  aspectHeatmap(forshow, scale='none', col=redgreen(64), main=main,
                ColSideColors=ColSideColors,
                RowSideColors=RowSideColors,
                Colv=Colv,
                Rowv=as.dendrogram(object@gc),
                ...
                )
  invisible(object)
})

if (!isGeneric("makeFigures"))
  setGeneric("makeFigures",
             function(object, DIR=NULL, ...) { standardGeneric("makeFigures") }
             )

setMethod("makeFigures", "Thresher", function(object, DIR=NULL, ...) {
  if (is.null(DIR)) opar <- par(ask=TRUE)
  fname <- gsub("\\.", "-", object@name) # LaTeX-safe
# fig1
  if(!is.null(DIR)) png(file.path(DIR, paste(fname, "-01-scree.png", sep="")),
               width=640, height=640, bg="white")
  screeplot(object)
  if(!is.null(DIR)) dev.off()
# fig2
  if(!is.null(DIR)) png(file.path(DIR, paste(fname, "-02-bayes.png", sep="")),
               width=800, height=600, bg="white")
  plot(object@ag, list(object@agfun), lwd=3, main=paste(object@name, "(Bayesian Anaysis)"))
  if(!is.null(DIR)) dev.off()
# fig3
  if(!is.null(DIR)) png(file.path(DIR, paste(fname, "-03-spca.png", sep="")),
               width=640, height=640, bg="white")
  scatter(object)
  if(!is.null(DIR)) dev.off()
# fig4 (a, b, ...)
  M <- max(1, trunc((1 + object@pcdim)/2))
  for (i in 1:M) {
    if(!is.null(DIR)) png(file.path(DIR,
                                    paste(fname, "-04", LETTERS[i],
                                          "-loadings.png", sep="")),
                          width=640, height=640, bg="white")
    plot(object, ij=c(2*i-1, 2*i))
    if(!is.null(DIR)) dev.off()
  }
# fig5
  if(!is.null(DIR)) png(file.path(DIR, paste(fname, "-05-heatmap.png", sep="")),
               width=1920, height=640, pointsize=12, bg="white")
  heat(object, wExp=3, margins=c(5,15))
  if(!is.null(DIR)) dev.off()
# finish
  if (is.null(DIR)) par(opar)
  invisible(object)
})

