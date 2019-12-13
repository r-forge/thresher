wd <<- choose.dir()
setwd(wd)

save.image(file='Idiograph1.rda')
load("Idiograph1.RData")
euroFreq <- OSU_AML

### from European/302-summaryFigs
### euroFreq must have a column named Chromosome
### format is chr1, chr2,... ect
makeIdiogram <- function(colname, color, res=1368) {
  par(mfrow=c(2,24), mai=c(0, 0.1, 1, 0.1), bg='white')
  for (I in c(1:22, "X", "Y")) {
    chrname <- paste('chr', I, sep='')
    dumbposn <- seq(1, 250000000, length=2500)
    clap <- cytobandLocations[cytobandLocations$Chromosome == chrname,]
    segset <- euroFreq[euroFreq$Chromosome == chrname,]
    y <- rep(NA, length(dumbposn))
    v <- rep(NA, length(dumbposn))
    for(i in 1:nrow(clap)) {
      y[clap[i, "loc.start"] <= dumbposn & 
          dumbposn <= clap[i, "loc.end"] ] <- as.numeric(clap[i, "Stain"])
      v[clap[i, "loc.start"] <= dumbposn & 
          dumbposn <= clap[i, "loc.end"] ] <- as.numeric(segset[i, colname])
    }
    par(mai=c(0, 0.1, 0.5, 0.05))
    image(1:1, dumbposn, matrix(rev(y), nrow=1), col=idiocolors, bty='n',
          xlab='', ylab='', xaxt='n', yaxt='n', zlim=c(1, 8), 
          main=paste("Chr", I))
    pts <- max(dumbposn)-c(min(clap$loc.start), max(clap$loc.end))
    abline(h=pts)
    lines(c(1.4, 1.4), pts)
    lines(c(0.6, 0.6), pts)
    par(mai=c(0, 0.01, 0.5, 0.05))
    barplot(rev(v), horiz=T, border=NA, col="blue",
            xlim=c(0, res/100), yaxs="i", xaxt="n", space=0)
  }
}

### from02-wicell-biIdiogram
### DATA must have three columns, one of which is "Chromosome" containing
### entries of the form "chr1".
biIdiogram  <- function(DATA, leftcol, rightcol, 
                        pal = c("blue", "red"), 
                        horizontal=FALSE, nrows = 2) {
  if(!nrows %in% 1:4) {
    stop("Number of rows must be 1, 2, 3, or 4.")
  }
  opar <- par(bg="white")
  on.exit(par(opar))
  # vertical layout, two rows
  L1 <- function() {
    layout(matrix(1:72, 1, 72, byrow = TRUE),
           heights = c(1), widths=rep(c(1.1, 0.8, 1.1), times = 24))
  }
  L2 <- function() {
    layout(matrix(1:72, 2, 36, byrow = TRUE),
           heights = c(1,1), widths=rep(c(1.1, 0.8, 1.1), times = 12))
  }
  L3 <- function() {
    layout(matrix(1:72, 3, 24, byrow = TRUE),
           heights = c(1,1, 1), widths=rep(c(1.1, 0.8, 1.1), times = 8))
  }
  L4 <- function() {
    layout(matrix(1:72, 4, 18, byrow = TRUE),
           heights = c(1,1, 1, 1), widths=rep(c(1.1, 0.8, 1.1), times = 6))
  }
  switch(nrows, L1(), L2(), L3(), L4())
  
  for (I in c(1:22, "X", "Y")) { # for each chromosome
    chrname <- paste('chr', I, sep='')
    dumbposn <- seq(1, 250000000, length=2500)
    clap <- cytobandLocations[cytobandLocations$Chromosome == chrname,]
    segset <- DATA[DATA$Chromosome == chrname,]
    leftvals <- rightvals <- y <- rep(NA, length(dumbposn))
    for(J in 1:nrow(clap)) {
      y[clap[J, "loc.start"] <= dumbposn & 
          dumbposn <= clap[J, "loc.end"] ] <- as.numeric(clap[J, "Stain"])
      rightvals[clap[J, "loc.start"] <= dumbposn & 
                  dumbposn <= clap[J, "loc.end"] ] <- as.numeric(segset[J, rightcol])
      leftvals[clap[J, "loc.start"] <= dumbposn & 
                 dumbposn <= clap[J, "loc.end"] ] <- as.numeric(segset[J, leftcol])
    }
    resn <- max(max(DATA[,leftcol]), max(segset[,rightcol]))
    ## left bars
    par(mai=c(0, 0.001, 0.5, 0.001))
    barplot(-rev(leftvals), horiz=T, border=NA, col=pal[1],
            xlim=c(-1.05*resn, 0), yaxs="i", xaxt="n", space=0)
    ## chromosomes
    par(mai=c(0, 0.02, 0.5, 0.02))
    image(1:1, dumbposn, matrix(rev(y), nrow=1), col=idiocolors, bty='n',
          xlab='', ylab='', xaxt='n', yaxt='n', zlim=c(1, 8),
          main=paste("Chr", I), cex=0.8)
    pts <- max(dumbposn) - c(min(clap$loc.start), max(clap$loc.end))
    abline(h=pts)
    lines(c(1.4, 1.4), pts)
    lines(c(0.6, 0.6), pts)
    ## right bars
    par(mai=c(0, 0.001, 0.5, 0.001))
    barplot(rev(rightvals), horiz=T, border=NA, col=pal[2],
            xlim=c(0, 1.05*resn), yaxs="i", xaxt="n", space=0)
  }
  layout(1,1,1)
  par(new = TRUE)
  legend(0.4, 400, c(leftcol, rightcol), col=pal, pch=15)
  #  box()
  invisible(DATA)
}


getFreq <- function(name, standard = wisc) {
  tenv <- new.env()
  on.exit(rm(tenv))
  fname <- paste(name, "LGF.Rda", sep='')
  load(file.path(lymPaths$clean, fname), envir = tenv)
  DATA <- get(paste(name, "Data", sep=""), envir = tenv)
  keep <- apply(DATA, 1, function(x) any(x > 0))
  DATA <- DATA[keep,]
  FREQ <- 100*apply(DATA, 2, mean)
  temp <- data.frame(Loss = FREQ[grep("loss", names(FREQ))],
                     Gain = FREQ[grep("gain", names(FREQ))],
                     Fusion = FREQ[grep("fusion", names(FREQ))])
  rownames(temp) <- sub("loss.", "", rownames(temp))
  w <- which(rownames(temp) == "19q13.2.")
  if (length(w) == 1) rownames(temp)[w] <- "19q13.2"
  temp[rownames(standard),]
}

# DATA must include a column called "Chromosome" containing
# entries of the form "chr1".
horizontalStacks  <- function(DATA, columns,
                              pal = "blue",   #mypal
                              nrows = 2) {
  if(!nrows %in% 1:4) {
    stop("Number of rows must be 1, 2, 3, or 4.")
  }
  N <- length(columns)
  if (N < 1) stop("You need to supply the name of at least one data column.")
  if (N > 10) stop("Unable to show more than ten data columns.")
  
  opar <- par(bg="white")
  on.exit(par(opar))
  # vertical layout, two rows
  corder <- function(N1) {
    M <- matrix(1:(24*N1), ncol=N1, byrow=TRUE)
    M <- M[, seq(ncol(M), 1, -1)]
    as.vector(t(M))
  }
  L1 <- function(ncolumn) {
    N1 <- ncolumn + 1
    layout(matrix(corder(N1), 24*N1, 1, byrow = FALSE),
           widths = c(1), heights=rep(c(0.8, rep(1.1, ncolumn)), times = 24))
  }
  L2 <- function(ncolumn) {
    N1 <- ncolumn + 1
    layout(matrix(corder(N1), 12*N1, 2, byrow = FALSE),
           widths = c(1,1), heights=rep(c(0.8, rep(1.1, ncolumn)), times = 12))
  }
  L3 <- function(ncolumn) {
    N1 <- ncolumn + 1
    layout(matrix(corder(N1), 8*N1, 3, byrow = FALSE),
           widths = c(1,1, 1), heights=rep(c(0.8, rep(1.1, ncolumn)), times = 8))
  }
  L4 <- function(ncolumn) {
    N1 <- ncolumn + 1
    layout(matrix(corder(N1), 6*N1, 4, byrow = FALSE),
           widths = c(1,1, 1, 1), heights=rep(c(0.8, rep(1.1, ncolumn)), times = 6))
  }
  
  switch(nrows, L1(N), L2(N), L3(N), L4(N))
  for (I in c(1:22, "X", "Y")) { # for each chromosome
    chrname <- paste('chr', I, sep='')
    dumbposn <- seq(1, 250000000, length=2500)
    clap <- cytobandLocations[cytobandLocations$Chromosome == chrname,]
    segset <- DATA[DATA$Chromosome == chrname,]
    resn <- max(max(DATA[,columns]))
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
    mtext(paste("Chr", I), side=2, at=1, las=2, line=0.5)
    pts <- c(min(clap$loc.start), max(clap$loc.end))
    abline(v=pts)
    lines(pts, c(1.4, 1.4))
    lines(pts, c(0.6, 0.6))
    ## right bars
    for (K in 1:length(columns)) {
      cname <- columns[K]
      vals <- NA*y
      for(J in 1:nrow(clap)) {
        vals[clap[J, "loc.start"] <= dumbposn & 
               dumbposn <= clap[J, "loc.end"] ] <- as.numeric(segset[J, cname])
      }
      par(mai=c(0.001, 0.6, 0.001, 0))
      barplot(vals, horiz=F, border=NA, col=pal[K],
              ylim=c(0, 1.05*resn), xaxs="i", yaxt="n", space=0)
    }
  }
  layout(1,1,1)
  par(new = TRUE)
  legend(0.78*par("usr")[2], 0.45*par("usr")[4], columns, col=pal, pch=15)
  invisible(DATA)
}

# DATA must include a column called "Chromosome" containing
# entries of the form "chr1".
stackIdiogram  <- function(DATA, columns,
                           pal = "blue",         #mypal,
                           horizontal=FALSE, nrows = 2) {
  if (horizontal) {
    return(horizontalStacks(DATA, coplumns, pal, nrows))
  }
  if(!nrows %in% 1:4) {
    stop("Number of rows must be 1, 2, 3, or 4.")
  }
  opar <- par(bg="white")
  on.exit(par(opar))
  # vertical layout, two rows
  L1 <- function(ncolumn) {
    N1 <- ncolumn + 1
    layout(matrix(1:(24*N1), 1, 24*N1, byrow = TRUE),
           heights = c(1), widths=rep(c(0.8, rep(1.1, ncolumn)), times = 24))
  }
  L2 <- function(ncolumn) {
    N1 <- ncolumn + 1
    layout(matrix(1:(24*N1), 2, 12*N1, byrow = TRUE),
           heights = c(1,1), widths=rep(c(0.8, rep(1.1, ncolumn)), times = 12))
  }
  L3 <- function(ncolumn) {
    N1 <- ncolumn + 1
    layout(matrix(1:(24*N1), 3, 8*N1, byrow = TRUE),
           heights = c(1,1, 1), widths=rep(c(0.8, rep(1.1, ncolumn)), times = 8))
  }
  L4 <- function(ncolumn) {
    N1 <- ncolumn + 1
    layout(matrix(1:(24*N1), 4, 6*N1, byrow = TRUE),
           heights = c(1,1, 1, 1), widths=rep(c(0.8, rep(1.1, ncolumn)), times = 6))
  }
  N <- length(columns)
  if (N < 1) stop("You need to supply the name of at least one data column.")
  if (N > 10) stop("Unable to show more than ten data columns.")
  switch(nrows, L1(N), L2(N), L3(N), L4(N))
  
  for (I in c(1:22, "X", "Y")) { # for each chromosome
    chrname <- paste('chr', I, sep='')
    dumbposn <- seq(1, 250000000, length=2500)
    clap <- cytobandLocations[cytobandLocations$Chromosome == chrname,]
    segset <- DATA[DATA$Chromosome == chrname,]
    resn <- max(max(DATA[,columns]))
    y <- rep(NA, length(dumbposn))
    for(J in 1:nrow(clap)) {
      y[clap[J, "loc.start"] <= dumbposn & 
          dumbposn <= clap[J, "loc.end"] ] <- as.numeric(clap[J, "Stain"])
    }
    ## chromosomes
    par(mai=c(0, 0.05, 0.6, 0.02))
    image(1:1, dumbposn, matrix(rev(y), nrow=1), col=idiocolors, bty='n',
          xlab='', ylab='', xaxt='n', yaxt='n', zlim=c(1, 8),
          cex=0.8)
    mtext(paste("Chr", I), side=3, at=1, las=2, line=0.5, cex=1)
    pts <- max(dumbposn) - c(min(clap$loc.start), max(clap$loc.end))
    abline(h=pts)
    lines(c(1.4, 1.4), pts)
    lines(c(0.6, 0.6), pts)
    ## right bars
    for (K in 1:length(columns)) {
      cname <- columns[K]
      vals <- NA*y
      for(J in 1:nrow(clap)) {
        vals[clap[J, "loc.start"] <= dumbposn & 
               dumbposn <= clap[J, "loc.end"] ] <- as.numeric(segset[J, cname])
      }
      par(mai=c(0, 0.001, 0.6, 0.001))
      barplot(rev(vals), horiz=T, border=NA, col=pal[K],
              xlim=c(0, 1.05*resn), yaxs="i", xaxt="n", space=0)      
    }
  }
  layout(1,1,1)
  par(new = TRUE)
  legend(0.3*par("usr")[2], 0.2*par("usr")[4], columns, col=pal, pch=15)
  invisible(DATA)
}

makeIdiogram()
biIdiogram(OSU_AML,3,5)
#graphs AML_LOSS on column 3, 6 gives it color
stackIdiogram(OSU_AML_LOSS,3,6) 
stackIdiogram(OSU_AML_GAIN,3)
stackIdiogram(OSU_AML_FUSION,3,3)
horizontalStacks(OSU_AML_GAIN,3)

