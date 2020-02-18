library(RCytoGPS)

set.seed(53106)
CL <- RCytoGPS:::cytobandLocations
fake <- as.data.frame(matrix(runif(nrow(CL) * 10), ncol = 10))
CL <- cbind(CL, fake) # extra colnames V1, ..., V10

plot(as.numeric(CL$Chromosome)) # check chromosome order
showGenome()

genomeBarplot(CL$V1, "brown")

opar <- par(mfrow = c(3,1))
genomeBarplot(CL$V2, "green")
genomeBarplot(CL$V3, "orange")
genomeBarplot(CL$V5, "purple")
par(opar)
