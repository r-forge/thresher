library(RCytoGPS)

set.seed(53106)
CL <- RCytoGPS:::cytobandLocations
fake <- as.data.frame(matrix(runif(nrow(CL) * 10), ncol = 10))
CL <- cbind(CL, fake) # extra colnames V1, ..., V10
all(diff(as.numeric(CL$Chromosome)) >= 0) # check chromosome order

CD <- CytobandData(CL)
barplot(CD, what = "V10")
barplot(CD, what = list("V1")) # automagically converts list to vector

try( barplot(CD, what = c("V1", "V10")) )    # fail
try( barplot(CD, what = list("V1", "V10")) ) # fail

# confirm usability within preset par grid
opar <- par(mfrow = c(3,1))
barplot(CD, "V2", "green")
barplot(CD, "V3", "orange")
barplot(CD, "V5", "purple")
par(opar)
