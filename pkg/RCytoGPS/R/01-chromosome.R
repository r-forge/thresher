setClass("Chromosome",
         slots = c(
           name = "character",
           label = "character",
           grid = "numeric",
           range = "numeric",
           stain = "numeric")
)

Chromosome <- function(I, res = 2500, maxbase = 250000000) {
  if (length(I) > 1) stop("Can only handle one chromosome at a time.")
  if (!(I %in% c(1:22, "X", "Y"))) stop("Invalid chromosome name!")
  chrname <- paste('chr', I, sep='')
  label <- paste('Chr', I)
  dumbposn <- seq(1, maxbase, length=res)
  YY <- cytobandLocations[cytobandLocations$Chromosome == chrname,]
  y <- rep(NA, length(dumbposn))
  for(i in 1:nrow(YY)) {
    y[YY[i, "loc.start"] <= dumbposn & 
      dumbposn <= YY[i, "loc.end"] ] <- as.numeric(YY[i, "Stain"])
  }
  new("Chromosome",
      name = chrname,
      label = label,
      grid = dumbposn,
      range = c(min(YY$loc.start), max(YY$loc.end)),
      stain = y)
}

setMethod("image", signature = "Chromosome", function(x,
            horiz = FALSE, mai = NULL, showBandNames = FALSE, ...) {
  if (is.null(mai)) {
    if (horiz) {
      mai <- c(0, 0.1, 0.5, 0.05)
    } else {
      mai <- c(0.05, 0.6, 0.02, 0)
    }
  }
  opar <- par(mai=mai)
  on.exit(par(opar))
  if (horiz) {
    pts <- max(x@grid) - x@range
    image(1:1, x@grid, matrix(rev(x@stain), nrow=1), col=idiocolors, bty='n',
          xlab='', ylab='', xaxt='n', yaxt='n', zlim=c(1, 8),
          main=x@label)
    rect(0.6, pts[1], 1.4, pts[2])
    if (showBandNames) {
      foo <- findLabels(x@name)
      text(1, max(x@grid) - foo$midpt, rownames(foo), col=foo$labCol)
    }
  } else {
    pts <- x@range
    image(x@grid, 1:1, matrix(x@stain, ncol=1), col=idiocolors, bty='n',
          xlab='', ylab='', xaxt='n', yaxt='n', zlim=c(1, 8),
          cex=0.8)
    mtext(x@label, side=2, at=1, las=2, line=0.5)
    rect(pts[1], 0.6, pts[2], 1.4)
    if (showBandNames) {
      foo <- findLabels(x@name)
      text(foo$midpt, 1, rownames(foo), col=foo$labCol, srt=90)
    }
  }
  invisible(x)
})

idioLabelColors <- c(acen = "white", gneg = "black", gpos100 = "white", gpos25 = "black",
                     gpos50 = "black", gpos75 = "white", gvar = "white", stalk = "black")

findLabels <- function(chrom) {
  if (length(chrom) > 1) {
    warning("More than two chromosomes supplied; only using the first one")
    chrom <- chrom[1]
  }
  if (!grepl("^chr", chrom)) { # add the prefix
    chrom <- paste("chr", chrom, sep = "")
  }
  X <- cytobandLocations[cytobandLocations$Chromosome == chrom,]
  X$midpt <- (X$loc.start + X$loc.end)/2
  X$labCol <- idioLabelColors[X$Stain]
  X
}
