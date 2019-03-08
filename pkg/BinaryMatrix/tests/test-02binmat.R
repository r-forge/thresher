library("BinaryMatrix")
data("CML500")
data(lgfFeatures)

### need a plain old binary matrix
working <- t(CML500@binmat)
### and a data.fame describing the feature-columns.
feat <- lgfFeatures[colnames(working),]

# construction
egg <- BinaryMatrix(working, feat)
summary(egg)

# remove extra identical feature vectors
egg <- removeDuplicateFeatures(egg)
summary(egg)
all( dim(egg) == c(511, 446) )

