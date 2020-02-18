library(RCytoGPS)

#windows(width=8, height=18)

set.seed(53106)
CL <- RCytoGPS:::cytobandLocations
fake <- as.data.frame(matrix(runif(nrow(CL) * 10), ncol = 10))
CL <- cbind(CL, fake) # extra colnames V1, ..., V10

plot1Chrom(CL, "V1", 1)
plot1Chrom(CL, c("V1", "V2"), 11)
plot1Chrom(CL, c("V1", "V2", "V3","V4", "V5", "V6"), 13)
plot1Chrom(CL, paste("V", 1:10, sep=""), 6)

opar <- par(mfrow=c(2,1))
plot1Chrom(CL, c("V1", "V2", "V3"), 2)
plot1Chrom(CL, c("V4", "V5", "V6"), 13)
par(opar)

colnames(CL)[6:14]  <- paste(rep(LETTERS[1:3], each=3),
                             rep(c("Loss", "Gain", "Fusion"), times = 3),
                             sep=".")
singles(CL, LETTERS[1:3], 3)
