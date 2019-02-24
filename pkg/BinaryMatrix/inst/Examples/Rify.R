source("parseBands.R")# to get L:GF annotations

library(rjson)

slurp <- readLines("CML_Mitelman_2018.json")
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

ID <- paste("S",
            unlist(lapply(wow, function(X) X[["id"]])),
            sep ='')

msg <- unlist(lapply(wow, function(X) X[["message"]]))
unique(msg) # so we don't need this

clone <- unlist(lapply(wow, function(X) X[["clone"]]))
indic <- unlist(lapply(wow, function(X) X[["indicator"]]))

fullID <- paste(ID, clone, indic, sep='.')
CMLInfo <- data.frame(ID = ID, Clone = clone, Indicator = indic)
rownames(CMLInfo) <- fullID
rm(ID, msg, clone, indic, fullID)

CMLData <- matrix(NA, nrow = length(wow), ncol = length(wow[[1]]$LGF))
dim(CMLData)
for (I in 1:length(wow)) {
  X <- wow[[I]]$LGF
  CMLData[I,] <- X
}
rm(I, X, wow)

(ncol(CMLData) == nrow(lgfFeatures))
colnames(CMLData) <- rownames(lgfFeatures)
rownames(CMLData) <- rownames(CMLInfo)
CMLData[1:10, 1:4]


save(CMLData, CMLInfo, lgfFeatures, file = "../../data/CML.rda")
