library("Mercator")

## test names part
set.seed(853845)
X <- matrix(sample(c(0,1), 10*30, replace=TRUE), 10, 30)
B <- BinaryMatrix(X)

J <- Mercator(B, "jaccard", "mds", 3)
J <- addVisualization(J, "tsne", perplexity=5)
J <- addVisualization(J, "graph")
G <- J@view[["graph"]]
## There used to be a bug that caused thie next line to fail when row names
## or column names were NULL. If this works, then we haven't broken it again.
plot(G$graph, layout=G$layouts[["nicely"]])
plot(J, view = "graph", layout = "nicely")
## There is an issue in igraph if nodes are placed at (nearly?) identical
## locations. This was fixed inside the 'plot' method  by 'jittering" the
## locations.
plot(J, view = "graph", layout = "mds")
plot(J, view = "graph", layout = "tsne")

## Need to see if we still assume that visualizations exist
J <- Mercator(B, "jaccard", "hclust", 3)
J <- addVisualization(J, "graph")
G <- J@view[["graph"]]
names(G$layouts)

