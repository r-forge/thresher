library("BinaryMatrix")
data("CML1000")

# transposition
bm <-t(CML1000)
bm <- removeDuplicateFeatures(bm)
summary(bm)

# run the thresher algorithm to remove "useless" vectors
thrash <- threshLGF(bm, cutoff = 0.5)
summary(thrash)

# convert back to base class
oval <- as(thrash, "BinaryMatrix")
slotNames(oval)
summary(oval)

