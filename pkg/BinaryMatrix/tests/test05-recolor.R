library("BinaryMatrix")
data("CML")
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
vis2 <- DistanceVis(oval, "sokal", "mds", K=8)
vis3 <- recolor(vis1, vis2)

#opar <- par(mfrow=c(1,3),  cex=2)
plot(vis1@view[[1]], col=vis1@colv, pch=vis1@symv)
plot(vis2@view[[1]], col=vis2@colv, pch=vis2@symv)
plot(vis3@view[[1]], col=vis3@colv, pch=vis3@symv)
#par(opar)
