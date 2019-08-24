library("Mercator")
data("CML500")

# Jaccard
vis1 <- Mercator(CML500, "jacc", "mds", K=20)
head(vis1@symv) # make sure names propagate correctly
head(vis1@colv)
plot(vis1) # default MDS view

barplot(vis1)
scatter(vis1, view = "mds")

vis1 <- addVisualization(vis1, "tsne", perplexity=30)
plot(vis1, view = "tsne")

scatter(vis1, view = "tsne", colramp = grDevices::terrain.colors)

vis1 <- addVisualization(vis1, "hclust")
plot(vis1, view = "hclust")

# now test the igraph part
vis1 <- addVisualization(vis1, "graph")
plot(vis1, view = "graph") # default layout
plot(vis1, view = "graph", layout = "nicely") # default layout
plot(vis1, view = "graph", layout = "mds")
plot(vis1, view = "graph", layout = "tsne")
