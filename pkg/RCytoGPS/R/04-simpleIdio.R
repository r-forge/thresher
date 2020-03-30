plot2Chrom <- function(DATA, leftcol, rightcol, chr,
                       pal = c("blue", "red"),
                       horiz = FALSE, axes = TRUE, legend = FALSE, resn = NULL,
                       sigcolumn = NA, sigcut = 0.01, alpha = 63) {
  ## check sigcolumn, sigcut, alpha
  if (!is.na(sigcolumn)) {
    if (length(sigcut) == 0) stop("You must supply at least one significance cutoff!")
    if (length (sigcut) != length(alpha)) stop("Lengths of 'sigcut' and 'alpha' do not match.")
    sigcut <- sort(sigcut)
    alpha <- sort(alpha)
  }
  ## check valid short chromosome name
  if ( !(chr %in% c(1:22, "X", "Y")) ) stop("Invalid chromosome number.")
  chrname <- paste("chr", chr, sep="")
  ## check valid stack height
  while(length(pal) < 2) pal <- c(pal, pal)
  ## check that all columns exist
  if ( !(leftcol %in% colnames(DATA)) )
    stop("Unrecognized (left) column name.")
  if ( !(rightcol %in% colnames(DATA)) )
    stop("Unrecognized (right) column name.")
  opar <- par(c("bg", "new", "mai", "usr"))
  on.exit(par(opar))
  par(bg = "white", mai=c(0,0,0,0))
  ## fake plot to white screen so we can use par(new=TRUE) in loop below
  plot(0, 0, xaxt="n", yaxt="n", xlab="", ylab="", type="n", axes=FALSE)
  ## get the figure size in inches
  fin  <- par("fin")
  ## set the resolution
  V0 <- 15
  V1 <- 45

  dumbposn <- seq(1, 250000000, length=2500)
  CL <- cytobandLocations
  clap <- CL[CL$Chromosome == chrname,]
  segset <- DATA[DATA$Chromosome == chrname,]
  if (is.null(resn)) {
    resn <- max(max(segset[, c(leftcol, rightcol)]))
  }
  y <- left <- right <- dex <- rep(NA, length(dumbposn))
  for(J in 1:nrow(clap)) {
    y[clap[J, "loc.start"] <= dumbposn &
      dumbposn <= clap[J, "loc.end"] ] <- as.numeric(clap[J, "Stain"])
    left[clap[J, "loc.start"] <= dumbposn &
         dumbposn <= clap[J, "loc.end"] ] <- as.numeric(segset[J, leftcol])
    right[clap[J, "loc.start"] <= dumbposn &
          dumbposn <= clap[J, "loc.end"] ] <- as.numeric(segset[J, rightcol])
    if (!is.na(sigcolumn)) {
      dex[clap[J, "loc.start"] <= dumbposn &
          dumbposn <= clap[J, "loc.end"] ] <- 1 + sum(segset[J, sigcolumn] < sigcut)
    } else {
      dex[clap[J, "loc.start"] <= dumbposn &
             dumbposn <= clap[J, "loc.end"] ] <- 1 + length(sigcut)
    }
  }
  leftpal <- c(makeTransparent(pal[1], alpha), pal[1])
  rightpal <- c(makeTransparent(pal[2], alpha), pal[2])
  if(horiz) {
    vres <- fin[2]/100
    hres <- fin[1]/(V0 + 2*V1)

    ## right
    par(new = TRUE,
        mai=c(vres, (1 + V0 + V1)*hres, 15*vres, 2*hres))
    barplot(rev(right), horiz=TRUE, border=NA, col = rev(rightpal[dex]),
            xlim=c(0, 1.05*resn), yaxs="i",
            space=0, axes = FALSE)
    if (axes) {
      axis(3) # on top
      U <- par("usr")
      W <- (U[1] + U[2])/2
      mtext(rightcol, at = W, side=3, line=2.5)
    }
    ## chromosome in the middle
    par(new=TRUE,
        mai=c(vres, (V1 + 2)*hres, 15*vres, (V1 + 2)*hres))
    image(Chromosome(chr), mai=par("mai"), horiz = TRUE)
    ## left, pointing backwards
    par(new = TRUE,
        mai=c(vres, 2*hres, 15*vres, (1 + V0 + V1)*hres))
    barplot(-rev(left), horiz=TRUE, border=NA, col = rev(leftpal[dex]),
            xlim=c(-1.05*resn, 0), yaxs="i",
            space=0, axes=FALSE)
    if (axes) {
      axis(3) # on top
      U <- par("usr")
      W <- (U[1] + U[2])/2
      mtext(leftcol, at = W, side=3, line=2.5)
    }
  } else {
    vres <- fin[2]/(V0 + 2*V1)
    hres <- fin[1]/100

    ## right goes on top, pointing up
    par(new = TRUE,
        mai=c((1 + V0 + V1)*vres, 10*hres, 2*vres, hres))
    barplot(right, border=NA, col = rightpal[dex],
            ylim=c(0, 1.05*resn), xaxs="i", ylab = rightcol,
            space=0, axes = axes)
    ## chromosome in the middle
    par(new=TRUE,
        mai=c((V1 + 2)*vres, 10*hres, (V1 + 2)*vres, hres))
    image(Chromosome(chr), mai=par("mai"), horiz = horiz)
    ## left goes on the bottom, and points down
    par(new = TRUE,
        mai=c(2*vres, 10*hres, (1 + V0 + V1)*vres, hres))
    barplot(-(left), border=NA, col = leftpal[dex],
            ylim=c(-1.05*resn, 0), xaxs="i", ylab = leftcol,
            space=0, axes = axes)
  }
  par(opar)
  invisible(DATA)
}

biIdiogram  <- function(DATA, leftcol, rightcol,
                        pal = c("blue", "red"), nrows = 2,
                        horiz = FALSE, axes = TRUE, legend = FALSE,
                        sigcolumn = NA, sigcut = 0.01, alpha = 63) {
  if(!nrows %in% 1:4) {
    stop("Number of rows must be 1, 2, 3, or 4.")
  }
  opar <- par(c("mfrow", "mai", "usr"))
  on.exit(par(opar))

  if (horiz) { # horizontal layout of vertical plots
    L1 <- function() { par(mfrow = c(24, 1)) }
    L2 <- function() { par(mfrow = c(12, 2)) }
    L3 <- function() { par(mfrow = c(8, 3)) }
    L4 <- function() { par(mfrow = c(6, 4)) }
  } else { # vertical layout of horizontal plots
    L1 <- function() { par(mfrow = c(1, 24)) }
    L2 <- function() { par(mfrow = c(2, 12)) }
    L3 <- function() { par(mfrow = c(3, 8)) }
    L4 <- function() { par(mfrow = c(4, 6)) }
  }
  switch(nrows, L1(), L2(), L3(), L4())

  resn <- max(max(DATA[, c(leftcol, rightcol)]))
  
  for (I in c(1:22, "X", "Y")) { # for each chromosome
    plot2Chrom(DATA, leftcol, rightcol, I, pal, !horiz,
               axes = axes, resn = resn,
               sigcolumn = sigcolumn, sigcut = sigcut, alpha = alpha)
  }
  if (legend) {
    par(opar)
    par(mai=c(1, 1, 1, 1), usr = c(0,1, 0, 1))
    if (horiz) {
      legend("right", c(leftcol, rightcol), col = pal, lwd = 5)
    } else {
      legend("bottom", c(leftcol, rightcol), col = pal, lwd = 5)
    }
  }
}
