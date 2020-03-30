library(RCytoGPS)
## simulate data
set.seed(53106)
CL <- RCytoGPS:::cytobandLocations
fake <- as.data.frame(matrix(runif(nrow(CL) * 10), ncol = 10))
fake$PV <- runif(nrow(fake), 0, 1)
CL <- cbind(CL, fake) # extra colnames V1, ..., V10
CD <- CytobandData(CL)

try( image(CD))     # fails, since you need 'chr'
try( image(CD, 3) ) # fails, since you also need 'what'

### windows(width=16, height=16)
## (plot1Chrom) one chromsome, stacks, vertical
image(CD, chr = 1, what = "V1", debug = TRUE)
image(CD, chr = 1, what = c("V1", "V2"), debug = TRUE)
image(CD, chr = 13, what = c("V1", "V2", "V3","V4", "V5", "V6"), debug = TRUE)
image(CD, chr = 6, what = paste("V", 1:10, sep=""), debug = TRUE)

## (plot1Chrom) one chromosome, stacks, horizontal
image(CD, 6, paste("V", 1:2, sep=""), horiz = TRUE, debug = TRUE)
image(CD, 6, paste("V", 1:2, sep=""), axes = FALSE, debug = TRUE)
image(CD, 6, paste("V", 1:2, sep=""), horiz = TRUE, axes = FALSE, debug = TRUE)

image(CD, 8, paste("V", 1:2, sep=""), horiz = FALSE,
      sigcolumn = "PV", sigcut = c(0.01, 0.05), alpha = c(127, 192))

image(CD, 8, paste("V", 1:2, sep=""), horiz = TRUE,
      sigcolumn = "PV", sigcut = c(0.01, 0.05), alpha = c(127, 192), debug = TRUE)

# combined (plot1Chrom)
opar <- par(mfrow=c(2,1))
image(CD, 2, c("V1", "V2", "V3"), debug = TRUE)
image(CD, 13, c("V4", "V5", "V6"), debug = TRUE)
par(opar)

## (makeIdiogram)
image(CD, "all", "V10", pal = "forestgreen", debug = TRUE) # no names  is better
image(CD, "all", "V10", pal = "forestgreen", debug = TRUE, axes = TRUE) # names
image(CD, "all", "V10", pal = "forestgreen", legend = TRUE, debug = TRUE)
image(CD, "all", "V10", pal = "forestgreen", debug = TRUE,
      sigcolumn = "PV", sigcut = c(0.01, 0.05), alpha = c(63, 155))

## (plot2Chrom) one chromosome, two-sided
image(CD, 3, list("V1", "V2"), debug = TRUE) # default axes = TRUE and horiz = FALSE
image(CD, 3, list("V1", "V2"), axes = FALSE, debug = FALSE)
image(CD, 3, list("V1", "V2"), horiz = TRUE, debug = TRUE)
image(CD, 3, list("V1", "V2"), horiz = TRUE, axes = FALSE, debug = TRUE) # still names?

image(CD, 8, list("V1", "V2"),
      sigcolumn = "PV", sigcut = c(0.01, 0.05), alpha = c(127, 192),
      debug = TRUE)
image(CD, 8, list("V1", "V2"), horiz = TRUE,
      sigcolumn = "PV", sigcut = c(0.01, 0.05), alpha = c(127, 192),
      debug = TRUE)

## (biIdiogram)
image(CD, chr = "all", list("V1", "V2"), horiz=TRUE, debug = TRUE) # default axes = FALSE
image(CD, chr = "all", list("V1", "V2"), horiz=FALSE, debug  = TRUE)
image(CD, chr = "all", list("V1", "V2"), horiz=TRUE, axes = TRUE, debug = TRUE)
image(CD, chr = "all", list("V1", "V2"), horiz=FALSE, axes = TRUE, debug  = TRUE)
image(CD, chr = "all", list("V1", "V2"), horiz=TRUE, nrows = 3)
image(CD, chr = "all", list("V1", "V2"), horiz=FALSE, nrows = 3)
image(CD, chr = "all", list("V1", "V2"), horiz=TRUE, nrows = 4)
image(CD, chr = "all", list("V1", "V2"), horiz=FALSE, nrows = 4)
image(CD, chr = "all", list("V1", "V2"), horiz=TRUE, nrows = 1)
image(CD, chr = "all", list("V1", "V2"), horiz=FALSE, nrows = 1)
## still need to test legend and sigcolumn
image(CD, chr = "all", list("V1", "V2"), horiz=TRUE, legend = TRUE, debug = TRUE)
image(CD, chr = "all", list("V1", "V2"), horiz=FALSE, legend = TRUE, debug  = TRUE)

image(CD, chr = "all", list("V1", "V2"),
      sigcolumn = "PV", sigcut = c(0.01, 0.05), alpha = c(127, 192),
      debug = TRUE)
image(CD, chr = "all", list("V1", "V2"), horiz = TRUE,
      sigcolumn = "PV", sigcut = c(0.01, 0.05), alpha = c(127, 192),
      debug = TRUE)

# stackIdiograms
mycolumns <- paste("V", 1:3, sep="")
image(CD, chr = "all", mycolumns, debug = TRUE)
image(CD, chr = "all", mycolumns, nrows = 3, debug = TRUE)
image(CD, chr = "all", mycolumns, nrows = 4, horiz=FALSE, debug = TRUE)
image(CD, chr = "all", mycolumns, nrows = 4, horiz=TRUE, debug = TRUE)
image(CD, chr = "all", mycolumns, nrows = 4, axes=TRUE, horiz=FALSE, debug = TRUE)
image(CD, chr = "all", mycolumns, nrows = 4, axes=TRUE, horiz=TRUE, debug = TRUE)
## still need to check axes, legend, and sigcolumn

image(CD, chr = "all", mycolumns, debug = TRUE, legend = TRUE) # no good place to put it
image(CD, chr = "all", mycolumns, horiz = TRUE, debug = TRUE, legend = TRUE) # no good place to put it


image(CD, chr = "all", mycolumns, horiz = TRUE,
      sigcolumn = "PV", sigcut = c(0.01, 0.05), alpha = c(63, 155),
      debug = TRUE)
