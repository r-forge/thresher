library("BinaryMatrix")
data("CML1000")
data(lgfFeatures)

bm <- BinaryMatrix(t(CML1000@binmat),
                   lgfFeatures[rownames(CML1000@binmat),])
bm <- removeDuplicateFeatures(bm)
summary(bm)

# jacc, pear, manh, euc
vis1 <- DistanceVis(bm, "jacc", "mds", K=20)
plot(vis1@view[[1]], col=vis1@colv, pch=vis1@symv)

vis1 <- addVisualization(vis1, "tsne", perplexity=20)
plot(vis1@view[[2]]$Y, col=vis1@colv, pch=vis1@symv)

vis1 <- addVisualization(vis1, "hclust")
plot(vis1@view[[3]])
