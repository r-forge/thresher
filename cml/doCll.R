ok <- system("perl jsonify.pl CLL_Mitelman_2018.txt")

source("parseBands.R")# to get LGF annotations

library(rjson)

slurp <- readLines("CLL_Mitelman_2018.json")
longline <- paste(slurp, collapse="")
wow <- fromJSON(longline)
length(wow)
wow$producer
wow$version
wow$source

wow <- wow$results
length(wow)
names(wow[[1]])

rm(slurp, longline)
gc()

ID <- paste("CLL",
            unlist(lapply(wow, function(X) X[["id"]])),
            sep ='')
msg <- unlist(lapply(wow, function(X) X[["message"]]))
unique(msg) # so we don't need this

clone <- unlist(lapply(wow, function(X) X[["clone"]]))
indic <- unlist(lapply(wow, function(X) X[["indicator"]]))

fullID <- paste(ID, clone, indic, sep='.')
any(duplicated(fullID))

CLLInfo <- data.frame(ID = ID, Clone = clone, Indicator = indic)
rownames(CLLInfo) <- fullID
rm(ID, msg, clone, indic, fullID)

CLLData <- matrix(NA, nrow = length(wow), ncol = length(wow[[1]]$LGF))
dim(CLLData)
for (I in 1:length(wow)) {
  X <- wow[[I]]$LGF
  CLLData[I,] <- X
}
rm(I, X, wow)

any(is.na(CLLData))
(ncol(CLLData) == nrow(lgfFeatures))
colnames(CLLData) <- rownames(lgfFeatures)
rownames(CLLData) <- rownames(CLLInfo)
CLLData[1:10, 1:4]


save(CLLData, CLLInfo, lgfFeatures, file = "CLL.rda")

# cleanup
file.remove("CLL_Mitelman_2018.json")
file.remove("lgfFeatures.rda")
