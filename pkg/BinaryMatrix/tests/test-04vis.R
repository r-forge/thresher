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
# convert back to base class
oval <- as(thrash, "BinaryMatrix")

# jacc, pear, manh, euc
vis1 <- DistanceVis(oval, "jacc", "mds", K=8)
plot(vis1@view[[1]], col=vis1@colv, pch=vis1@symv)

vis1 <- addVisualization(vis1, "tsne", perplexity=10)
plot(vis1@view[[2]]$Y, col=vis1@colv, pch=vis1@symv)

vis1 <- addVisualization(vis1, "hclust")
plot(vis1@view[[3]])
