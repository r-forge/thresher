library(RCytoGPS)

try( x <- Chromosome(0) )   # fail
try( x <- Chromosome("Z") ) # fail
try( x <- Chromosome(1:2) ) # fail
try (x <- Chromosome(c(1, "Z"))) # fail
try (x <- Chromosome(list(1, "Z"))) # fail

x <- Chromosome(11)
image(x)
image(x, showBandNames = TRUE)
image(x, horiz = TRUE)
image(x, horiz = TRUE, showBandNames = TRUE)

image(x, col="red") # silently ignores bad input?
