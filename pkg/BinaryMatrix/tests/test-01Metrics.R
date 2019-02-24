library("BinaryMatrix")
data("CML")

# load("../data/CML.Rda")
set.seed(92547)
working <- CMLData[sample(nrow(CMLData), 1100),
                   sample(ncol(CMLData), 276)]
zero <- apply(working, 1, function(x) all(x == 0))
table(zero)
working <- working[!zero,]
zero <- apply(working, 2, function(x) all(x == 0))
table(zero)
working <- working[, !zero]
dim(working)

jacc <- binaryDistance(working, "jaccard")
all(dim(as.matrix(jacc)) == c(250, 250))
plot(hclust(jacc, "ward.D2"))

canb <- binaryDistance(working, "canberra")
all(dim(as.matrix(canb)) == c(250, 250))
plot(hclust(canb, "ward.D2"))

X <- binaryDistance(working, "sokal")
all(dim(as.matrix(X)) == c(250, 250))
plot(hclust(X, "ward.D2"))

X <- binaryDistance(working, "hamming")
all(dim(as.matrix(X)) == c(250, 250))
plot(hclust(X, "ward.D2"))

X <- binaryDistance(working, "russell")
all(dim(as.matrix(X)) == c(250, 250))
plot(hclust(X, "ward.D2"))

X <- binaryDistance(working, "goodman")
all(dim(as.matrix(X)) == c(250, 250))
plot(hclust(X, "ward.D2"))

X <- binaryDistance(working, "manhattan")
all(dim(as.matrix(X)) == c(250, 250))
plot(hclust(X, "ward.D2"))

X <- binaryDistance(working, "binary")
all(dim(as.matrix(X)) == c(250, 250))
plot(hclust(X, "ward.D2"))

X <- binaryDistance(working, "euclid")
all(dim(as.matrix(X)) == c(250, 250))
plot(hclust(X, "ward.D2"))

