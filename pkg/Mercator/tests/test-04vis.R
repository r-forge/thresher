library("Mercator")
data("CML500")

## Jaccard
vis1 <- Mercator(CML500, "jacc", "mds", K=20)
head(Mercator:::symv(vis1)) # make sure names propagate correctly
head(Mercator:::colv(vis1))
head(vis1@clusters)
plot(vis1, main = "MDS") # default MDS view

barplot(vis1)
scatter(vis1, view = "mds", main = "MDS")

## t-SNE
vis1 <- addVisualization(vis1, "tsne", perplexity=30)
plot(vis1, view = "tsne", main = "t-SNE")

scatter(vis1, view = "tsne", colramp = grDevices::terrain.colors, main = "t-SNE")

## hclust
vis1 <- addVisualization(vis1, "hclust")
plot(vis1, view = "hclust")

## now test the igraph part
vis1 <- addVisualization(vis1, "graph")
plot(vis1, view = "graph") # default layout
plot(vis1, view = "graph", layout = "nicely") # default layout
plot(vis1, view = "graph", layout = "mds")
plot(vis1, view = "graph", layout = "tsne")

## umap
vis1 <- addVisualization(vis1, "umap")
plot(vis1, view = "umap", main="UMAP")

scatter(vis1, view = "umap", main="UMAP")

## som
vis1 <- addVisualization(vis1, "som")
plot(vis1, view = "som")
plot(vis1, view = "som", type = "counts")
plot(vis1, view="som", type = "dist.neighbours", palette.name = topo.colors)
