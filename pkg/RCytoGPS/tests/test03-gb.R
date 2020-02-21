library(RCytoGPS)

set.seed(53106)
CL <- RCytoGPS:::cytobandLocations
fake <- as.data.frame(matrix(runif(nrow(CL) * 10), ncol = 10))
CL <- cbind(CL, fake) # extra colnames V1, ..., V10

plot(as.numeric(CL$Chromosome)) # check chromosome order
showGenome()

CD <- CytobandData(CL)
barplot(CD, what = "V10")

opar <- par(mfrow = c(3,1))
barplot(CD, "V2", "green")
barplot(CD, "V3", "orange")
barplot(CD, "V5", "purple")
par(opar)
