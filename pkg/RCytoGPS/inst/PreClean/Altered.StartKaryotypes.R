library(tidyverse)
library(stringr)
x <- read.csv("startKaryotypes.csv") # original data
y <- read.delim("badStrings.txt") # list of abnormalities strings need to delete from z

z <- data.frame(x)#copy of original data
View(y)#double check what I need to get rid of

#relacing the strings of the abnormality
z[,4:5] <- lapply(z[,4:5], function(x) gsub(".ish","",x, fixed = TRUE))
z[,4:5] <- lapply(z[,4:5], function(x) gsub("/nonclonal w","",x, fixed = TRUE))
z[,4:5] <- lapply(z[,4:5], function(x) gsub("/nonclonal[1]","",x, fixed = TRUE))
z[,4:5] <- lapply(z[,4:5], function(x) gsub("/nonclonal w/clonal abnormalities[2]","",x, fixed = TRUE))
z[,4:5] <- lapply(z[,4:5], function(x) gsub("[3, one is 4n]","",x, fixed = TRUE))
z[,4:5] <- lapply(z[,4:5], function(x) gsub("/nonclonal[2]","",x, fixed = TRUE))
z[,4:5] <- lapply(z[,4:5], function(x) gsub("clonal","",x, fixed = TRUE))
z[,4:5] <- lapply(z[,4:5], function(x) gsub("/4n[1]","",x, fixed = TRUE))
z[,4:5] <- lapply(z[,4:5], function(x) gsub("/4n[2]","",x, fixed = TRUE))
z[,4:5] <- lapply(z[,4:5], function(x) gsub("/4n[3]","",x, fixed = TRUE))
z[,4:5] <- lapply(z[,4:5], function(x) gsub("/4n[4]","",x, fixed = TRUE))
z[,4:5] <- lapply(z[,4:5], function(x) gsub("(D7S486+;D7Z1+,D7S486-)","",x, fixed = TRUE))
z[,4:5] <- lapply(z[,4:5], function(x) gsub("46,XX[2]","",x, fixed = TRUE))
z[,4:5] <- lapply(z[,4:5], function(x) gsub("(D13S319-,LAMP+)","",x, fixed = TRUE))
z[,4:5] <- lapply(z[,4:5], function(x) gsub("/46,XX[17]","",x, fixed = TRUE))

View(x)#viewing original data to be compared
View(z)#Viewing altered data to see if replacing worked

write.csv(z, file = "Altered.StartKaryotypes.csv")
