library(Thresher)

# unstructured matrix
set.seed(9948489)
dumb <- matrix(rnorm(100*12), ncol=12)
colnames(dumb) <- paste("G", 1:12, sep='')
thresh <- Thresher(dumb)
# test 'failure' version of constructor
r <- new("Reaper", thresh,
         useLoadings=FALSE,
         keep=FALSE,
         nGroups=NA,
         fit=NA,
         bic=NA,
         allfits=list(),
         metric='pearson')
# now construct it for real
reap <- Reaper(thresh)
screeplot(reap)
plot(reap)
scatter(reap)
heat(reap)

# matrices with structure
set.seed(250264)
rho <- 0.5
nProtein <- 16
splinter <- sample((nProtein/2) + (-3:3), 1)
sigma1 <- matrix(rho, ncol=nProtein, nrow=nProtein)
diag(sigma1) <- 1
sigma2 <- sigma1
sigma2[(1+splinter):nProtein, 1:splinter] <- 0
sigma2[1:splinter, (1+splinter):nProtein] <- 0
nonsense <- matrix(rnorm(nProtein^2, 0, 0.001), nProtein)
sigma2 <- sigma2 + nonsense + t(nonsense)
# now simulate the data
thresh <- SimThresher(sigma2, nSample=300)
summary(thresh@delta)
# create Reaper
reap <- Reaper(thresh)
reap@pcdim # two real components
screeplot(reap, col='gold', lcol='black')
reap@nGroups # should be just two feature-clusters
plot(reap)
scatter(reap)
heat(reap)

colsch <- Thresher:::.makeColorScheme(4)
bin.hc  <- reap@signalSet@binaryClusters
bin.csc <- colsch[cutree(bin.hc, k=4)]
con.hc  <- reap@signalSet@continuousClusters
con.csc <- colsch[cutree(con.hc, k=4)]
heat(reap, Colv=as.dendrogram(bin.hc), ColSideColors=con.csc,
       main=paste(reap@name, "binary signals, continuous colors"))
heat(reap, Colv=as.dendrogram(con.hc), ColSideColors=bin.csc,
       main=paste(reap@name, "continuous signals, binary colors"))

if(FALSE) {
  makeFigures(reap)
}
