library(RCytoGPS)

sf <- system.file("Examples/JSONfiles", package = "RCytoGPS")
dir(sf)

## Read one explicitly named JSON file
temp <- readLGF("CytoGPS_Result1.json", folder = sf)
class(temp)
names(temp)
temp$source
class(temp$raw)
length(temp$raw)
names(temp$raw)
temp$raw[[1]]$Status
summary(temp$frequency)
temp$size

## Read all (well, two) JSON files fom a folder
temp <- readLGF(folder = sf)
class(temp)
names(temp)
temp$source
class(temp$raw)
length(temp$raw)
names(temp$raw)
temp$raw[[1]]$Status
temp$raw[[2]]$Status
summary(temp$frequency)
temp$size

## move about
home  <- getwd()
setwd(sf)
temp <- readLGF( "CytoGPS_Result1.json")
temp <- readLGF( "CytoGPS_Result1.json", folder = ".")
setwd("..")
temp <- readLGF(folder = "JSONfiles")
setwd(home)
