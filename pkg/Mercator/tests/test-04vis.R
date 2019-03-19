library("Mercator")
data("CML500")

# Jaccard
vis1 <- Mercator(CML500, "jacc", "mds", K=20)
head(vis1@symv) # make sure names propagate correctly
head(vis1@colv)
plot(vis1@view[[1]], col=vis1@colv, pch=vis1@symv)

vis1 <- addVisualization(vis1, "tsne", perplexity=30)
plot(vis1@view[[2]]$Y, col=vis1@colv, pch=vis1@symv)

vis1 <- addVisualization(vis1, "hclust")
plot(vis1@view[[3]])

# now test the igraph part
vis1 <- addVisualization(vis1, "graph")
G <- vis1@view[["graph"]]
plot(G$graph, layout=G$layouts[["nicely"]])
plot(G$graph, layout=G$layouts[["mds"]])
plot(G$graph, layout=G$layouts[["tsne"]])
