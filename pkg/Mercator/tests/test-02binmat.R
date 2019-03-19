library("Mercator")
data("CML500")

### need a plain old binary matrix
working <- t(CML500@binmat)
feat <- CML500@columnInfo
sams <- CML500@rowInfo
dim(working)
dim(feat)
dim(sams)

plain <- BinaryMatrix(working) # default row and column info
colonly <- BinaryMatrix(working, sams)

# construction, using both info parts
egg <- BinaryMatrix(working, sams, feat)
summary(egg)

# remove extra identical feature vectors
egg <- removeDuplicateFeatures(egg)
summary(egg)
all( dim(egg) == c(511, 446) )

# transposing back using the t-operator.
gge <- t(egg)
summary(gge)
