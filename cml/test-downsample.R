library("BinaryMatrix")
sessionInfo()

### get the CML data set, and remove duplicate LGF features
load("CML.rda")
dumb <- data.frame(Name = rownames(CMLData))
rownames(dumb) <- rownames(CMLData)
bam <- new("BinaryMatrix",
           binmat = CMLData,
           columnInfo = lgfFeatures,
           rowInfo = CMLInfo,
           info = list(),
           history =  "Manual")

bm <- BinaryMatrix(CMLData, lgfFeatures, CMLInfo)
bm <- removeDuplicateFeatures(bm)

### create the transposed binary matrix (columns = samples)
flipped <- t(bm)
rm(CMLData, CMLInfo, lgfFeatures, bm)

### Compute two visualizations based on Jaccard distance
### This takes a fair amount oif time, since there are more
### than 5000 samples to cluster.
vis <- DistanceVis(flipped, "jacc", "mds", K=48)
plot(vis@view[[1]], col=vis@colv, pch=vis@symv)
### Add a t-SNE plot, which is also slow
vis <- addVisualization(vis, "tsne")
plot(vis@view[[2]]$Y, col=vis@colv, pch=vis@symv, cex=1.5)

### Test the "downsample" routine. This requires a distance matrix
### (like Jaccard) along with a parameter to define which neighbors
### are "close". For the latter, we are going to use the 15th
### percentile of distances.
mn <- apply(DM <- as.matrix(vis@distance), 1, function(x) min(x[x>0]))
summary(mn)
R <- quantile(mn, c(0.10, 0.15, 0.20))
# prep
set.seed(58667)
pick1000 <- downsample(1000, DM, R[2])
summary(pick1000)

### First, we use the (new) subset operator on the DistanceVis
### class to get the smaller set of samples with the "inherited"
### visualizations. We didn't create a hierarchical clustering
### dendrogram above because subsetting is also really super-slow.
wow <- vis[pick1000]
plot(wow@view[[1]], col=wow@colv, pch=wow@symv, cex=1.5)
plot(wow@view[[2]]$Y, col=wow@colv, pch=wow@symv, cex=1.5)
### Notice that the *center* of the large magenta "eye" is absent.


### Now we use the subset operator on the BinaryMatrix class
### to create a smaller matrix.
CML1000 <- flipped[, pick1000]
### And we create new visualizations
vis1000 <- DistanceVis(CML1000, "jacc", "mds", K=20)
vis1000 <- addVisualization(vis1000, "tsne")

windows(width=14, height=7)
par(mfrow=c(1,2))
plot(vis1000@view[[1]], col=vis1000@colv, pch=vis1000@symv,
     cex=1.5, main="MDS, New Colors")
plot(vis1000@view[[1]], col=vis@colv[pick1000], pch=vis@symv[pick1000],
     cex =1.5, main="MDS, Old Colors")

plot(vis1000@view[[2]]$Y, col=vis1000@colv, pch=vis1000@symv,
     cex=1.5, main="TSNE, New Colors")
plot(vis1000@view[[2]]$Y, col=vis@colv[pick1000], pch=vis@symv[pick1000],
     cex=1.5, main="TSNE, Old Colors")

### For use in other examples, we downsample further
set.seed(235564)
pick500 <- downsample(500, DM, R[2])
CML500 <- flipped[, pick500]


save(CML1000, file = "CML1000.rda")
save(CML500, file = "CML500.rda")
