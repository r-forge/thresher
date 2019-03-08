library("BinaryMatrix")
data("CML500")

# remove extra identical feature vectors
CML500 <- removeDuplicateFeatures(CML500)
# jacc, pear, manh, euc
vis1 <- DistanceVis(CML500, "jacc", "mds", K=8)
vis2 <- DistanceVis(CML500, "sokal", "mds", K=8)
vis3 <- recolor(vis1, vis2)

# par(mfrow=c(1,3),  cex=1.5)
plot(vis1@view[[1]], col=vis1@colv, pch=vis1@symv)
plot(vis2@view[[1]], col=vis2@colv, pch=vis2@symv)
plot(vis3@view[[1]], col=vis3@colv, pch=vis3@symv)

A <- BinaryMatrix:::recoverCluster(vis2@colv, vis2@symv)
B <- BinaryMatrix:::recoverCluster(vis3@colv, vis3@symv)
table(A, B)
X <- BinaryMatrix:::recoverCluster(vis1@colv, vis1@symv)
table(A, X)
table(B, X)
