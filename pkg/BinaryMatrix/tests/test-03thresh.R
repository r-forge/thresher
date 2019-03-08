library("BinaryMatrix")
data("CML1000")
data(lgfFeatures)

bm <- BinaryMatrix(t(CML1000@binmat),
                   lgfFeatures[rownames(CML1000@binmat),])
bm <- removeDuplicateFeatures(bm)
summary(bm)

# run the thresher algorithm to remove "useless" vectors
thrash <- threshLGF(bm, cutoff = 0.5)
summary(thrash)

# convert back to base class
oval <- as(thrash, "BinaryMatrix")
slotNames(oval)
summary(oval)

