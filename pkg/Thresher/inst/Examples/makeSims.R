#############################################555555
### Used to simulate a set of related examples that
### can illustrate Thresher and Reaper methods.
library(Thresher)

set.seed(98576) # for reproducibility

### save eveything in lists
savedSims <- list() # SimThresher objects
savedReap <- list() # parallel Reaper objects

# parameters
nProtein <- 20
splinter <- sample((nProtein/2) + (-3:3), 1)
positive <- sample(nProtein, nProtein/2)
negative <- (1:nProtein)[!((1:nProtein) %in% positive)]
posi <- positive[positive <= splinter]
nega <- negative[negative <= splinter]
# one group, unsigned
rho <- rnorm(1, 0.5, 0.1)
sigma1 <- matrix(rho, ncol=nProtein, nrow=nProtein)
diag(sigma1) <- 1
# two groups, unsigned
sigma2 <- sigma1
sigma2[(1+splinter):nProtein, 1:splinter] <- 0
sigma2[1:splinter, (1+splinter):nProtein] <- 0
# one group, signed
sigma3 <- sigma1
sigma3[positive, negative] <- -sigma3[positive, negative]
sigma3[negative, positive] <- -sigma3[negative, positive]
# two groups, signed
sigma4 <- sigma2
sigma4[positive, negative] <- -sigma4[positive, negative]
sigma4[negative, positive] <- -sigma4[negative, positive]
# two groups, mixed
sigma5 <- sigma2
sigma5[posi, nega] <- -sigma5[posi, nega]
sigma5[nega, posi] <- -sigma5[nega, posi]
# number of samples
nSample <- round(rnorm(1, 300, 60))

### clean up
print(splinter)
print(rho)
rm(nega, negative, posi, positive, rho, splinter)

### Actually perfomr the simulations.
counter <- 0
for (nm in paste("sigma", 1:5, sep='')) {
# add noise
  ss <- matrix(0, nProtein+2, nProtein+2)
  diag(ss) <- 1
  ss[1:nProtein, 1:nProtein] <- get(nm)
# basic analysis, which runs Threshedr on simulated data
  counter <- counter+1
  value <- SimThresher(ss, nSample, 
                       paste(nm, counter, sep="."),
                       method='auer.gervini')
  savedSims[[counter]] <- value
# Now apply Reaper
  reaped <- Reaper(value, useLoadings=TRUE)
  savedReap[[counter]] <- reaped
}


if (FALSE) {
  for (i in 1:5)
    makeFigures(savedSims[[i]])

  for (i in 1:5)
    makeFigures(savedReap[[i]])
}

sigma <- list(sigma1, sigma2, sigma3, sigma4, sigma5)

save(savedSims, savedReap, sigma, file = "../../data/savedSims.Rda")

