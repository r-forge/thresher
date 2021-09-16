plot1Chrom <- function(DATA, columns,  chr, labels = columns,
                       pal = palette(), horiz = FALSE, axes=TRUE,
                       legend = FALSE, resn = NULL,
                       sigcolumn = NA, sigcut = 0.01, alpha = 63,
                       dlim = NULL) {
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
  NC <- length(columns)
  if (NC < 1) stop("You need to supply the name of at least one data column.")
  if (NC > 10) stop("Unable to show more than ten data columns.")
  while(length(pal) < NC) pal <- c(pal, pal)
  ## check that all columns exist
  if ( !all(columns %in% colnames(DATA)) ) stop("Unrecognized column name.")
  ## get the labels ready
  if (is.null(labels)) labels <- rep("", NC)
  while (length(labels) < NC) labels <- c(labels, labels)
  opar <- par(c("new", "mai", "bg", "usr"))
  on.exit(par(opar))
  ## fake plot to white screen so we can use par(new=TRUE) in loop below
  par(bg = "white", mai = c(0,0,0,0))
  plot(0, 0, xaxt="n", yaxt="n", xlab="", ylab="", type="n", axes=FALSE)
  ## get the figure size in inches
  fin <- par("fin")
  ## create a "vertical resolution" near 200
  V0     <- c(25, 25, 15, 10, 17, 14, 13, 12, 11, 10)
  V1     <- c(75, 50, 30, 20, 33, 28, 25, 22, 20, 18)
  if (horiz) {
    vres <- fin[2]/100
    hres <- fin[1]/(V0[NC] + NC*V1[NC])
  } else {
    vres <- fin[2]/(V0[NC] + NC*V1[NC])
    hres <- fin[1]/100
  }

  dumbposn <- seq(1, 250000000, length=2500)
  CL <- cytobandLocations
  clap <- CL[CL$Chromosome == chrname,]
  segset <- DATA[DATA$Chromosome == chrname,]
  if (is.null(resn)) {
    resn <- max(max(segset[, columns]))
  }
  y <- rep(NA, length(dumbposn))
  for(J in 1:nrow(clap)) {
    y[clap[J, "loc.start"] <= dumbposn & 
      dumbposn <= clap[J, "loc.end"] ] <- as.numeric(clap[J, "Stain"])
  }
  ## data bars
  for (II in 1:length(columns)) {
    K <- NC - II + 1
    column <- columns[K]
    kcolors <- c(makeTransparent(pal[K], alpha), pal[K])
    vals <- NA * y
    shades <- NA * y
    for(J in 1:nrow(clap)) {
      vals[clap[J, "loc.start"] <= dumbposn & 
           dumbposn <= clap[J, "loc.end"] ] <- as.numeric(segset[J, column])
      if (!is.na(sigcolumn)) {
        shades[clap[J, "loc.start"] <= dumbposn & 
               dumbposn <= clap[J, "loc.end"] ] <-
                    kcolors[1 + sum(segset[J, sigcolumn] < sigcut)]
      } else {
        shades[clap[J, "loc.start"] <= dumbposn & 
               dumbposn <= clap[J, "loc.end"] ] <- pal[K]
      }
    }
    if (horiz) {
      par(new = TRUE,
          mai=c(vres, hres * (1 + V0[NC] + (K-1)*V1[NC]),
                10*vres, hres * (1 + (II-1)*V1[NC])))
      barplot(rev(vals), horiz = horiz, border=NA, col = rev(shades), 
              xlim=c(0, 1.05*resn), yaxs="i", ylim = dlim,
              space=0, axes=FALSE)
      if (axes) {
        axis(3) # on top
        U <- par("usr")
        W <- (U[1] + U[2])/2
        mtext(labels[K], at = W, side=3, line=2.5)
      }
    } else {
      par(new = TRUE,
          mai=c(vres * (1 + V0[NC] + (K-1)*V1[NC]), 10*hres,
                vres * (1 + (II-1)*V1[NC]), hres))
      barplot(vals, horiz = horiz, border=NA, col = shades,
              ylim=c(0, 1.05*resn), xaxs="i", ylab = labels[K], xlim = dlim,
              space=0, axes = axes)
    }
  }
  ## chromosome
  if (horiz) {
    par(new=TRUE,
        mai=c(vres, 2*hres, 10*vres, hres * (1 + NC*V1[NC])))
  } else {
    par(new=TRUE,
        mai=c(2*vres, 10*hres, vres * (1 + NC*V1[NC]), hres))
  }
  image(Chromosome(chr), mai=par("mai"), horiz = horiz)
  par(opar)
  invisible(DATA)
}

makeIdiogram <- function(DATA, colname, color, axes = TRUE, legend = FALSE,
                       sigcolumn = NA, sigcut = 0.01, alpha = 63) {
  opar <- par(c("mfrow", "mai", "usr", "bg"))
  on.exit(par(opar))

  resn  <- max(max(DATA[, colname]))

  par(mfrow=c(2,12), mai=c(0, 0.1, 1, 0.1), bg='white')
  for (I in c(1:22, "X", "Y")) {
    plot1Chrom(DATA, colname, chr = I, pal = color,  horiz=TRUE, axes = axes,
               resn = resn,
               sigcolumn = sigcolumn, sigcut = sigcut, alpha = alpha)
  }
  if (legend) {
    par(opar)
    par(mai=c(1, 1, 1, 1), usr = c(0,1, 0, 1))
    box(lty="blank")
    legend("bottom", colname, col = color, lwd = 5)
  }
}

stackIdiogram  <- function(DATA, columns, pal = palette(), nrows = 2,
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

  resn  <- max(max(DATA[, columns]))

  for (I in c(1:22, "X", "Y")) { # for each chromosome
    plot1Chrom(DATA, columns, chr = I, pal = pal, horiz = !horiz, axes = axes,
               resn = resn,
               sigcolumn = sigcolumn, sigcut = sigcut, alpha = alpha)
  }
  if (legend) {
    par(opar)
    par(mai=c(1, 1, 1, 1), usr = c(0,1, 0, 1), new=TRUE)
    box(lty = "blank") # weird graphics hack
    if (horiz) {
      legend("right", columns, col = pal, lwd = 5)
    } else {
      legend("bottom", columns, col = pal, lwd = 5)
    }
  }
}
