library(RCytoGPS)

#windows(width=10, height=18)

set.seed(53106)
CL <- RCytoGPS:::cytobandLocations
fake <- as.data.frame(matrix(runif(nrow(CL) * 10), ncol = 10))
CL <- cbind(CL, fake) # extra colnames V1, ..., V10

plot1Chrom(CL, "V1", 3)
plot2Chrom(CL, "V1", "V2", 3)
plot2Chrom(CL, "V1", "V2", 3, horiz = TRUE)

biIdiogram(CL, "V1", "V2", horiz=TRUE)
biIdiogram(CL, "V1", "V2", horiz=FALSE)
biIdiogram(CL, "V1", "V2", horiz=TRUE, nrows = 3)
biIdiogram(CL, "V1", "V2", horiz=FALSE, nrows = 3)
biIdiogram(CL, "V1", "V2", horiz=TRUE, nrows = 4)
biIdiogram(CL, "V1", "V2", horiz=FALSE, nrows = 4)
biIdiogram(CL, "V1", "V2", horiz=TRUE, nrows = 1)
biIdiogram(CL, "V1", "V2", horiz=FALSE, nrows = 1)
