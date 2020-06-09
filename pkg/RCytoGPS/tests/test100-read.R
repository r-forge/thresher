library(RCytoGPS)

sf <- system.file("Examples/JSONfiles", package = "RCytoGPS")
dir(sf)

## Read one explicitly named JSON file
temp <- readLGF("CytoGPS_Result1.json", folder = sf)
class(temp)
length(temp)
temp <- temp[[1]]
temp$Status
dim(lgf  <-  temp$LGF)
lgf[, 1:4]
dot <- RCytoGPS:::Idioformat(lgf)
head(dot)

## Read all (well, two) JSON files fom a folder
temp <- readLGF(folder = sf)
length(temp)
temp[[1]]$Status
temp[[2]]$Status
temp[[1]]$LGF[,1:4]
temp[[2]]$LGF[,1:4]

