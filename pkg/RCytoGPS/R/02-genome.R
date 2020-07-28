showGenome <- function(altcol = "#FED4C4") {
  N<- as.numeric(cytobandLocations$Chromosome) %% 2
  M <- cumsum(table(cytobandLocations$Chromosome))
  midpt <- (c(M) + c(0, M[-length(M)]))/2
  labs <- sub("chr", "", names(midpt))

  image(1:length(N), 1:1, as.matrix(N), col=c("white", altcol),
        xaxt="n", xlab="", yaxt="n", ylab="")
  text(midpt, rep(1, length(midpt)), labs, adj=0.5)
}


genomeBarplot <- function(V, col = "blue", altcol = "#FED4C4", ylab="Percent", h = NULL) {
  ## get the figure size in inches
  fin = par("fin")
  ## define intenal "resolution" based on figure size
  HT0 <- 16
  HT1 <- 80
  vres <- fin[2]/(HT0 + HT1)
  hres <- fin[1]/1000
  opar <- par(c("new", "mai", "bg"))
  on.exit(par(opar))

#  par(bg = "white", mai=c(0.01, 0.5412, 0.1, 0.2772))
  par(bg = "white", mai=c(17*vres, 54*hres, 3*vres, 28*hres))
  barplot(V, border=NA, space=0, ylab=ylab, col=col, yaxs="i", xaxt="n")
  if (!is.null(h)) { abline(h = h) }
#  par(new=TRUE, mai=c(0.39, 0.86, 0.01, 0.63))
  par(new=TRUE, mai=c(10*vres, 86*hres, 81*vres, 63*hres))
  showGenome(altcol = altcol)
}
