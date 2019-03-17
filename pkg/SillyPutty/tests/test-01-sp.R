library(SillyPutty)

set.seed(3568)
nObj <- 110
fake <- matrix(rnorm(50*nObj), nrow=nObj, ncol = 50)
dimnames(fake) <- list(paste("Sample", 1:nObj, sep=''),
                       paste("Feature", 1:50, sep=''))
dis <- dist(fake)
labels <- sample(6, nObj, replace = TRUE)
names(labels) <- rownames(fake)
x <- SillyPutty(labels, dis)
plot(x@silhouette, col=1:6)

mds <- cmdscale(dis)
plot(mds, pch=16, col=x@cluster, cex=2)

set.seed(8732)
y <- RandomSillyPutty(dis, 6, N = 100)

sr <- summarizeRepeats(y)
