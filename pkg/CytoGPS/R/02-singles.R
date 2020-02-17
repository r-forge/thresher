plot1Chrom <- function(DATA, column,  chr) {
  if ( !(chr %in% c(1:22, "X", "Y")) ) stop("Invalid chromosome number.")
  if ( !(column %in% colnames(DATA)) ) stop("Unrecognized column name.")
  mai <- par("mai")
  fin = par("fin")
  mai[1] <- fin[2]/5
  opar <- par(bg = "white", mai = mai)
  on.exit(par(opar))

  dumbposn <- seq(1, 250000000, length=2500)
  clap <- cytobandLocations[cytobandLocations$Chromosome == chrname,]
  segset <- DATA[DATA$Chromosome == chrname,]
  resn <- max(max(segset[, column]))
  y <- rep(NA, length(dumbposn))
  for(J in 1:nrow(clap)) {
    y[clap[J, "loc.start"] <= dumbposn & 
      dumbposn <= clap[J, "loc.end"] ] <- as.numeric(clap[J, "Stain"])
  }
  ## chromosomes
  par(mai=c(0.05, 0.6, 0.02, 0))
  image(dumbposn, 1:1, matrix(y, ncol=1), col=idiocolors, bty='n',
        xlab='', ylab='', xaxt='n', yaxt='n', zlim=c(1, 8),
        cex=0.8)
  mtext(chrname, side=2, at=1, las=2, line=0.5)
  pts <- c(min(clap$loc.start), max(clap$loc.end))
  abline(v=pts)
  lines(pts, c(1.4, 1.4))
  lines(pts, c(0.6, 0.6))
  ## right bars
  for (K in 1:length(columns)) {
    cname <- column[K]
    vals <- NA*y
    for(J in 1:nrow(clap)) {
      vals[clap[J, "loc.start"] <= dumbposn & 
           dumbposn <= clap[J, "loc.end"] ] <- as.numeric(segset[J, cname])
    }
    par(mai=c(0.001, 0.6, 0.001, 0))
    barplot(vals, horiz=F, border=NA, col=pal[K],
            ylim=c(0, 1.05*resn), xaxs="i", ylab=paste("%", I),
            space=0)
  }
}



singles  <- function(DATA, columns, chr,pal = mypal) {
  N <- length(columns)
  if (N < 1) stop("You need to supply the name of at least one data column.")
  if (N > 10) stop("Unable to show more than ten data columns.")

  opar <- par(bg="white")
  on.exit(par(opar))
  # vertical layout
  corder <- function(N1) {
    M <- matrix(1:(3*N1), ncol=N1, byrow=TRUE)
    M <- M[, seq(ncol(M), 1, -1)]
    as.vector(t(M))
  }

  ncolumn <- length(columns)
  N1 <- ncolumn + 1
  # show in L-G-F order, one chromosome
  chrname <- paste('chr', chr, sep='')

  layout(matrix(corder(N1), 3*N1, 1, byrow = FALSE),
         widths = c(1), heights=rep(c(0.8, rep(1.1, ncolumn)), times = 3))
  for (I in c("Loss", "Gain", "Fusion")) { # for each abnormality type
    currcol <- paste(columns, I, sep=".")
    dumbposn <- seq(1, 250000000, length=2500)
    clap <- cytobandLocations[cytobandLocations$Chromosome == chrname,]
    segset <- DATA[DATA$Chromosome == chrname,]
    resn <- max(max(segset[,currcol]))
    y <- rep(NA, length(dumbposn))
    for(J in 1:nrow(clap)) {
      y[clap[J, "loc.start"] <= dumbposn & 
          dumbposn <= clap[J, "loc.end"] ] <- as.numeric(clap[J, "Stain"])
    }
    ## chromosomes
    par(mai=c(0.05, 0.6, 0.02, 0))
    image(dumbposn, 1:1, matrix(y, ncol=1), col=idiocolors, bty='n',
          xlab='', ylab='', xaxt='n', yaxt='n', zlim=c(1, 8),
          cex=0.8)
    mtext(chrname, side=2, at=1, las=2, line=0.5)
    pts <- c(min(clap$loc.start), max(clap$loc.end))
    abline(v=pts)
    lines(pts, c(1.4, 1.4))
    lines(pts, c(0.6, 0.6))
    ## right bars
    for (K in 1:length(columns)) {
      cname <- currcol[K]
      vals <- NA*y
      for(J in 1:nrow(clap)) {
        vals[clap[J, "loc.start"] <= dumbposn & 
          dumbposn <= clap[J, "loc.end"] ] <- as.numeric(segset[J, cname])
      }
      par(mai=c(0.001, 0.6, 0.001, 0))
      barplot(vals, horiz=F, border=NA, col=pal[K],
              ylim=c(0, 1.05*resn), xaxs="i", ylab=paste("%", I),
              space=0)
    }
  }
  layout(1,1,1)
  par(new = TRUE)
#  legend(0.78*par("usr")[2], 0.45*par("usr")[4], columns, col=pal, pch=15)
  invisible(DATA)
}

