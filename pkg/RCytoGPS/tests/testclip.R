library(RCytoGPS)
## simulate data
set.seed(53106)
CL <- RCytoGPS:::cytobandLocations
fake <- as.data.frame(matrix(runif(nrow(CL) * 10), ncol = 10))
fake$PV <- runif(nrow(fake), 0, 1)
CL <- cbind(CL, fake) # extra colnames V1, ..., V10
CD <- CytobandData(CL)

## clipping
image(CD, 7, list("V1", "V2"))
image(CD, 7, list("V1", "V2"), clip = TRUE)

RCytoGPS:::plot2Chrom(CD@DATA, "V1", "V2", chr = 7, horiz = TRUE, clip = TRUE)
RCytoGPS:::plot2Chrom(CD@DATA, "V1", "V2", chr = 7, horiz = FALSE, clip = TRUE)
image(Chromosome(7, maxbase = 1.6e8), horiz = TRUE)
