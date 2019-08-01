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
plot(G$graph, layout=G$layouts[["mds"]])
plot(G$graph, layout=G$layouts[["tsne"]])

## Need to see if we still assume that visualizations exist
J <- Mercator(B, "jaccard", "hclust", 3)
J <- addVisualization(J, "graph")
G <- J@view[["graph"]]
names(G$layouts)

