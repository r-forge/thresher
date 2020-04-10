library(RCytoGPS)

sf <- system.file("Examples/JSONfiles", package = "RCytoGPS")
dir(sf)

## Read one explicitly named JSON file
temp <- readLGF("CytoGPS_Result1.json", folder = sf)
class(temp)
length(temp)
temp <- temp[[1]]
dim(temp)
temp[, 1:4]

## Read all (well, two) JSON files fom a folder
temp <- readLGF(folder = sf)
length(temp)
lapply(temp, dim)
temp[[1]][,1:4]
temp[[2]][,1:4]
