setClass("Chromosome",
         slots = c(
           name = "character",
           label = "character",
           grid = "numeric",
           range = "numeric",
           stain = "numeric")
)

Chromosome <- function(I, res = 2500, maxbase = 250000000) {
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
            horizontal = FALSE, mai = NULL, ...) {
  if (is.null(mai)) {
    if (horizontal) {
      mai <- c(0.05, 0.6, 0.02, 0)
    } else {
      mai <- c(0, 0.1, 0.5, 0.05)
    }
  }
  opar <- par(mai=mai)
  on.exit(par(opar))
  if (horizontal) {
    pts <- x@range
    image(x@grid, 1:1, matrix(x@stain, ncol=1), col=idiocolors, bty='n',
          xlab='', ylab='', xaxt='n', yaxt='n', zlim=c(1, 8),
          cex=0.8)
    mtext(x@label, side=2, at=1, las=2, line=0.5)
    abline(v=pts)
    lines(pts, c(1.4, 1.4))
    lines(pts, c(0.6, 0.6))
  } else {
    pts <- max(x@grid) - x@range
    image(1:1, x@grid, matrix(rev(x@stain), nrow=1), col=idiocolors, bty='n',
          xlab='', ylab='', xaxt='n', yaxt='n', zlim=c(1, 8),
          main=x@label)
    abline(h=pts)
    lines(c(1.4, 1.4), pts)
    lines(c(0.6, 0.6), pts)
  }
})
