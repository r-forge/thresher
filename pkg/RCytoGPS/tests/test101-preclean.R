library(RCytoGPS)

sf <- system.file("PreClean", package = "RCytoGPS")
dir(sf)

bad <- read.delim(file.path(sf, "badStrings.txt"), header=FALSE)
bad <- as.vector(as.matrix(bad))
input <- read.csv(file.path(sf, "startKaryotypes.csv"))
colnames(input)
output <- read.csv(file.path(sf, "Altered.StartKaryotypes.csv"))[, -1]

myclean <- preclean(input, 4:5, bad)

summary(output[,4] == myclean[,4])
summary(output[,5] == myclean[,5])

if(FALSE) {
  foo <- data.frame(IN = input[,4], OUT = output[,4], MY = myclean[,4])[output[,4] != myclean[,4],]
  write.csv(foo, file="foo.csv")
}
