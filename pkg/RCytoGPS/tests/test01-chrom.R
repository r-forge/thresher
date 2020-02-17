library(RCytoGPS)

try( x <- Chromosome(0) )   # fail?
try( x <- Chromosome("Z") ) # fail

x <- Chromosome(11)
image(x)
image(x, horizontal = TRUE)
