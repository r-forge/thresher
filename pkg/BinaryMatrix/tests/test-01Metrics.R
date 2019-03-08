library("BinaryMatrix")
data("CML500")
temp <- removeDuplicateFeatures(CML500)
### need a plain old binary matrix
working <- temp@binmat
dim(working)
N <- dim(working)[2]

jacc <- binaryDistance(working, "jaccard")
all(dim(as.matrix(jacc)) == c(N, N))
plot(hclust(jacc, "ward.D2"))

canb <- binaryDistance(working, "canberra")
all(dim(as.matrix(canb)) == c(N, N))
plot(hclust(canb, "ward.D2"))

X <- binaryDistance(working, "sokal")
all(dim(as.matrix(X)) == c(N, N))
plot(hclust(X, "ward.D2"))

X <- binaryDistance(working, "hamming")
all(dim(as.matrix(X)) == c(N, N))
plot(hclust(X, "ward.D2"))

X <- binaryDistance(working, "russell")
all(dim(as.matrix(X)) == c(N, N))
plot(hclust(X, "ward.D2"))

X <- binaryDistance(working, "goodman")
all(dim(as.matrix(X)) == c(N, N))
plot(hclust(X, "ward.D2"))

X <- binaryDistance(working, "manhattan")
all(dim(as.matrix(X)) == c(N, N))
plot(hclust(X, "ward.D2"))

X <- binaryDistance(working, "binary")
all(dim(as.matrix(X)) == c(N, N))
plot(hclust(X, "ward.D2"))

X <- binaryDistance(working, "euclid")
all(dim(as.matrix(X)) == c(N, N))
plot(hclust(X, "ward.D2"))

