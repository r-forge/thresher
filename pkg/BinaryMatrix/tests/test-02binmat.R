library("BinaryMatrix")
data("CML")

# load("../data/CML.Rda")

# prep
set.seed(29574)
F <- sample(ncol(CMLData), 276)
working <- CMLData[sample(nrow(CMLData), 1100), F]

# construction
egg <- BinaryMatrix(working, lgfFeatures[F,])
summary(egg)

# remove extra identical feature vectors
egg <- removeDuplicateFeatures(egg)
summary(egg)
all( dim(egg) == c(1100, 159) )
