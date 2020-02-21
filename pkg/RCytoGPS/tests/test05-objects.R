library(RCytoGPS)

## windows(width=10, height=18)

set.seed(53106)
CL <- RCytoGPS:::cytobandLocations
fake <- as.data.frame(matrix(runif(nrow(CL) * 10), ncol = 10))
CL <- cbind(CL, fake) # extra colnames V1, ..., V10
CD <- CytobandData(CL)

## one chromsome, stacks, vertical
image(CD, what = "V1", chr = 1)
image(CD, what = c("V1", "V2"), chr = 1)
image(CD, what = c("V1", "V2", "V3","V4", "V5", "V6"), chr = 13)
image(CD, what = paste("V", 1:10, sep=""), chr = 6)

## windows(width=18, height=10)
# one chromsome, stacks, horizontal
plot1Chrom(CL, paste("V", 1:2, sep=""), 6, horiz = TRUE)

# combined
opar <- par(mfrow=c(2,1))
image(CD, c("V1", "V2", "V3"), 2)
image(CD, c("V4", "V5", "V6"), 13)
par(opar)

image(CD, "V10", chr = "all", pal = "forestgreen")

# one chromosome, two-sided
image(CD, list("V1", "V2"), 3)
image(CD, list("V1", "V2"), 3, horiz = TRUE)

# biIdiograms
image(CD, list("V1", "V2"), chr = "all", horiz=TRUE)
image(CD, list("V1", "V2"), chr = "all", horiz=FALSE)
image(CD, list("V1", "V2"), chr = "all", horiz=TRUE, nrows = 3)
image(CD, list("V1", "V2"), chr = "all", horiz=FALSE, nrows = 3)
image(CD, list("V1", "V2"), chr = "all", horiz=TRUE, nrows = 4)
image(CD, list("V1", "V2"), chr = "all", horiz=FALSE, nrows = 4)
image(CD, list("V1", "V2"), chr = "all", horiz=TRUE, nrows = 1)
image(CD, list("V1", "V2"), chr = "all", horiz=FALSE, nrows = 1)

# stackIdiograms
mycolumns <- paste("V", 1:6, sep="")
image(CD, mycolumns, chr = "all")
image(CD, mycolumns, chr = "all", nrows = 3)
image(CD, mycolumns, chr = "all", nrows = 4, horiz=FALSE)
image(CD, mycolumns, chr = "all", nrows = 4, horiz=TRUE)


## need to fix H-V
