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
r@nGroups
# now construct it for real
reap <- Reaper(thresh)
reap@nGroups
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
# now simulate the data
thresh <- SimThresher(sigma2, nSample=300)
thresh@nSample # 300
thresh@rho     # 0.5
thresh@pcdim   # 2
min(thresh@delta) > 0.5 # TRUE
# create Reaper
reap <- Reaper(thresh)
reap@pcdim   # two real components
reap@nGroups # and two clusters
screeplot(reap, col='gold', lcol='black')
plot(reap)
scatter(reap)
heat(reap)

colsch <- Thresher:::.makeColorScheme(4)
bin.hc  <- reap@signalSet@binaryClusters
bin.csc <- colsch[cutree(bin.hc, k=4)]
con.hc  <- reap@signalSet@continuousClusters
con.csc <- colsch[cutree(con.hc, k=4)]
heat(reap, Colv=as.dendrogram(bin.hc), ColSideColors=con.csc,
       main=paste(reap@name, "binary signals, cont. colors"))
heat(reap, Colv=as.dendrogram(con.hc), ColSideColors=bin.csc,
       main=paste(reap@name, "cont. signals, binary colors"))

if(FALSE) {
  makeFigures(reap)
}
