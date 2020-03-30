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

## windows(width=10, height=18)

## (plot1Chrom) one chromsome, stacks, vertical
image(CD, chr = 1, what = "V1")
image(CD, chr = 1, what = c("V1", "V2"))
image(CD, chr = 13, what = c("V1", "V2", "V3","V4", "V5", "V6"))
image(CD, chr = 6, what = paste("V", 1:10, sep=""))

## windows(width=18, height=10)
# (plot1Chrom) one chromosome, stacks, horizontal
image(CD, 6, paste("V", 1:2, sep=""), horiz = TRUE)

image(CD, 8, paste("V", 1:2, sep=""), horiz = TRUE,
      sigcolumn = "PV", sigcut = c(0.01, 0.05), alpha = c(75, 150))

image(CD, 8, paste("V", 1:2, sep=""), horiz = FALSE,
      sigcolumn = "PV", sigcut = c(0.01, 0.05), alpha = c(75, 150))

# combined
opar <- par(mfrow=c(2,1))
image(CD, 2, c("V1", "V2", "V3"))
image(CD, 13, c("V4", "V5", "V6"))
par(opar)

## (makeIdiogram)
image(CD, "all", "V10", pal = "forestgreen")
image(CD, "all", "V10", pal = "forestgreen", legend = TRUE)
image(CD, "all", "V10", pal = "forestgreen",
      sigcolumn = "PV", sigcut = c(0.01, 0.05), alpha = c(63, 155))


# one chromosome, two-sided
image(CD, 3, list("V1", "V2"))
image(CD, 3, list("V1", "V2"), horiz = TRUE)

# biIdiograms
image(CD, chr = "all", list("V1", "V2"), horiz=TRUE)
image(CD, chr = "all", list("V1", "V2"), horiz=FALSE)
image(CD, chr = "all", list("V1", "V2"), horiz=TRUE, nrows = 3)
image(CD, chr = "all", list("V1", "V2"), horiz=FALSE, nrows = 3)
image(CD, chr = "all", list("V1", "V2"), horiz=TRUE, nrows = 4)
image(CD, chr = "all", list("V1", "V2"), horiz=FALSE, nrows = 4)
image(CD, chr = "all", list("V1", "V2"), horiz=TRUE, nrows = 1)
image(CD, chr = "all", list("V1", "V2"), horiz=FALSE, nrows = 1)

# stackIdiograms
mycolumns <- paste("V", 1:6, sep="")
image(CD, chr = "all", mycolumns)
image(CD, chr = "all", mycolumns, nrows = 3)
image(CD, chr = "all", mycolumns, nrows = 4, horiz=FALSE)
image(CD, chr = "all", mycolumns, nrows = 4, horiz=TRUE)


## need to fix H-V
