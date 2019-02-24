library("BinaryMatrix")
data("CML")

# load("../data/CML.Rda")

# prep
set.seed(29574)
F <- sample(ncol(CMLData), 276)
working <- CMLData[sample(nrow(CMLData), 1100), F]

# construction
egg <- BinaryMatrix(working, lgfFeatures[F,])
# remove extra identical feature vectors
egg <- removeDuplicateFeatures(egg)

# run the thresher algorithm to remove "useless" vectors
thrash <- threshLGF(egg, cutoff = 0.5)
summary(thrash)

# convert back to base class
oval <- as(thrash, "BinaryMatrix")
slotNames(oval)
summary(oval)

