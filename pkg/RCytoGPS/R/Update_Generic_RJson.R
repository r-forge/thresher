---
  title: "Parsing Mitelman"
author: "Dwayne Tally"
date: "5/22/2019"
output:
  html_document:
  theme: yeti
highlight: kate
toc: true
editor_options: 
  chunk_output_type: inline
---
  
  # Get some data
  ```{r load.json}


#note in the vignette that turning the count into a percentage doesn't change the idiogram
#when downsampling, do NOT use a normal karyotype 

#possible solution, do something with pathing directory
#I.E. folder within the r package 
#generalize path structure, pull out any file at the end of the pathway
#localize the pathway from r package itself that makes the pathway and load in any json file from that directory
#call the dataframe the same name as the json file read in
#assume any json file that is within taht pathway
#give an error message that says the format isn't right
#make an empty file called input.jasonfile.location
#make a generic pathway location that can work on anymachine and use json file within that location
#cytoegps.org // to use json files from there
# Give the input file name to the function.

#setwd("C:/Users/firem/Documents/OSU/MyStuff")
# Load the package required to read JSON files.
if (!require("rjson")) {
  install.packages("rjson")
  library("rjson")
}
library(jsonlite)

wd <<- choose.dir()
setwd(wd)
message(sprintf("Current working dir: %s\n", wd))


filenames <- list.files(".", pattern="*.json", full.names=TRUE) # this should give you a character vector, with each file name represented by an entry

#doesn't work if you load in more than one json file
myJSON <- lapply(filenames, function(x) rjson::fromJSON(file=filenames)) # a list in which each element is one of your original JSON files

```

```{r output}
if (!require("splitstackshape")) {
  install.packages("splitstackshape")
  library("splitstackshape")
}
if (!require("dplyr")) {
  install.packages("dplyr")
  library("dplyr")
}
library(jsonlite)
Status <- sapply(myJSON[[1]]$output, "[", "status") #pulls out all the status object list

ID <- unlist(sapply(myJSON[[1]]$output, function(x){ #function to extract all the cell numbers/temporary ID
    lapply(x$parsing_result,"[[","cell_number")
}))


clone <- unlist(lapply(myJSON[[1]]$output, function(x){ #Function to parse all the loss data from parsing_result
  if((x$status)=="Success"){ 
    1:length(lapply(x$parsing_result, "[[", "loss_gain_fusion_computing"))
  }
}))

#ID <-  sapply(myJSON[[1]]$output[[1]], "[[", "cell_number")#pulls ID from second level

bands <- myJSON[[1]]$iscn2016_bands#pulls bands from top level

index <- sapply(myJSON[[1]]$output, function(x){ #function to get to the level of the karyotype index
  if((x$status)=="Success") x$parsing_result
})
cell_number <- lapply(j, function(x) length(x))
n.index <- lapply(index, function(x) length(x)) #function to essentially pull the karyotype/index from the data
n.index1 <- as.numeric(n.index)#converts it to numeric so we could manipulate the data


lgf <- lapply(myJSON[[1]]$output, function(x){ #Function to parse all the loss data from parsing_result
  if((x$status)=="Success"){ 
    lgf <- lapply(x$parsing_result, "[[","loss_gain_fusion_computing")
  }
})

loss1 <- sapply(lgf, function(x){
  lapply(x, "[[","loss")
})

gain1 <- sapply(lgf, function(x){
  lapply(x, "[[","gain")
})

fusion1 <- sapply(lgf, function(x){
  lapply(x, "[[","fusion")
})
id.index <- 1:length(lgf[[2]])
ID <- 2

loop, go into karyotype write everything to file and add one to sampling and then move on to the next one

#gain1 <- lapply(myJSON[[1]]$output, function(x){ #Function to parse all the gain data from parsing_result
 # if((x$status)=="Success") x$parsing_result[[1]]$loss_gain_fusion_computing$gain
#})

#fusion1 <- lapply(myJSON[[1]]$output, function(x){ #Function to parse all the fusion data from parsing_result
#  if((x$status)=="Success") x$parsing_result[[1]]$loss_gain_fusion_computing$fusion
#})

```


```{r seperate}
df_Loss <- as.data.frame(matrix(unlist(loss1), ncol = 916, byrow = TRUE), stringsAsFactors = FALSE)

colnames(df_Loss) <- c("Loss")
sep_Loss1<- cSplit(df_Loss,"Loss", sep = ",")#sepeartes by ,
colnames(sep_Loss1) <- bands#Gives the column names their associated bands
nw.df_Loss1 <- data.frame(sep_Loss1, check.names = FALSE, stringsAsFactors = FALSE) #combines the ID and sep_Loss1 ID= ID,
clean_Loss1 <- na.omit(nw.df_Loss1)#removes the NULL/NA rows along with their associated ID

#converts the dataframe to a matrix, to manipulate the data
m.clean_Loss1 <- as.matrix(clean_Loss1)
#converts c(0 or 0) to 0
m.clean_Loss1[m.clean_Loss1 =="c(0"] <- 0 
m.clean_Loss1[m.clean_Loss1 =="c(1"] <- 1 
m.clean_Loss1[m.clean_Loss1 =="0)"] <- 0 
#converts nonzeros to 1
m.clean_Loss1[m.clean_Loss1 =="1)"] <- 1 
m.clean_Loss1[m.clean_Loss1 ==2] <- 1 
m.clean_Loss1[m.clean_Loss1 =="c(2"] <- 1
m.clean_Loss1[m.clean_Loss1 =="2)"] <- 1
m.clean_Loss1[m.clean_Loss1 ==3] <- 1 
m.clean_Loss1[m.clean_Loss1 =="c(3"] <- 1
m.clean_Loss1[m.clean_Loss1 =="3)"] <- 1
m.clean_Loss1[m.clean_Loss1 ==4] <- 1 
m.clean_Loss1[m.clean_Loss1 =="c(4"] <- 1
m.clean_Loss1[m.clean_Loss1 =="4)"] <- 1
m.clean_Loss1[m.clean_Loss1 ==5] <- 1 
m.clean_Loss1[m.clean_Loss1 =="c(5"] <- 1
m.clean_Loss1[m.clean_Loss1 =="5)"] <- 1
m.clean_Loss1[m.clean_Loss1 ==6] <- 1 
m.clean_Loss1[m.clean_Loss1 =="c(6"] <- 1
m.clean_Loss1[m.clean_Loss1 =="6)"] <- 1
m.clean_Loss1[m.clean_Loss1 ==7] <- 1 
m.clean_Loss1[m.clean_Loss1 =="c(7"] <- 1
m.clean_Loss1[m.clean_Loss1 =="7)"] <- 1
m.clean_Loss1[m.clean_Loss1 ==8] <- 1 
m.clean_Loss1[m.clean_Loss1 =="c(8"] <- 1
m.clean_Loss1[m.clean_Loss1 =="8)"] <- 1
m.clean_Loss1[m.clean_Loss1 ==9] <- 1 
m.clean_Loss1[m.clean_Loss1 =="c(9"] <- 1
m.clean_Loss1[m.clean_Loss1 =="9)"] <- 1


df_Gain <- as.data.frame(matrix(unlist(gain1), ncol = 916, byrow = TRUE), stringsAsFactors = FALSE)
colnames(df_Gain) <- c("Gain")
sep_Gain1<- cSplit(df_Gain,"Gain", sep = ",")
colnames(sep_Gain1) <- bands
nw.df_Gain1 <- data.frame(sep_Gain1, check.names = FALSE, stringsAsFactors = FALSE)
clean_Gain1 <- na.omit(nw.df_Gain1)
#converts the dataframe to matrix
m.clean_Gain1 <- as.matrix(clean_Gain1)
#converts c(0 or 0) to 0
m.clean_Gain1[m.clean_Gain1 =="c(0"] <- 0
m.clean_Gain1[m.clean_Gain1 =="c(1"] <- 1
m.clean_Gain1[m.clean_Gain1 =="0)"] <- 0
#converts nonzeros to 1
m.clean_Gain1[m.clean_Gain1 =="1)"] <- 1
m.clean_Gain1[m.clean_Gain1 ==2] <- 1
m.clean_Gain1[m.clean_Gain1 =="c(2"] <- 1
m.clean_Gain1[m.clean_Gain1 =="2)"] <- 1
m.clean_Gain1[m.clean_Gain1 ==3] <- 1 
m.clean_Gain1[m.clean_Gain1 =="c(3"] <- 1
m.clean_Gain1[m.clean_Gain1 =="3)"] <- 1
m.clean_Gain1[m.clean_Gain1 ==4] <- 1 
m.clean_Gain1[m.clean_Gain1 =="c(4"] <- 1
m.clean_Gain1[m.clean_Gain1 =="4)"] <- 1
m.clean_Gain1[m.clean_Gain1 ==5] <- 1 
m.clean_Gain1[m.clean_Gain1 =="c(5"] <- 1
m.clean_Gain1[m.clean_Gain1 =="5)"] <- 1
m.clean_Gain1[m.clean_Gain1 ==6] <- 1 
m.clean_Gain1[m.clean_Gain1 =="c(6"] <- 1
m.clean_Gain1[m.clean_Gain1 =="6)"] <- 1
m.clean_Gain1[m.clean_Gain1 ==7] <- 1
m.clean_Gain1[m.clean_Gain1 =="c(7"] <- 1
m.clean_Gain1[m.clean_Gain1 =="7)"] <- 1
m.clean_Gain1[m.clean_Gain1 ==8] <- 1 
m.clean_Gain1[m.clean_Gain1 =="c(8"] <- 1
m.clean_Gain1[m.clean_Gain1 =="8)"] <- 1
m.clean_Gain1[m.clean_Gain1 ==9] <- 1 
m.clean_Gain1[m.clean_Gain1 =="c(9"] <- 1
m.clean_Gain1[m.clean_Gain1 =="9)"] <- 1

df_Fusion <- as.data.frame(matrix(unlist(fusion1), ncol = 916, byrow = TRUE), stringsAsFactors = FALSE)
colnames(df_Fusion) <- c("Fusion")
sep_Fusion1<- cSplit(df_Fusion,"Fusion", sep = ",")
colnames(sep_Fusion1) <- bands
nw.df_Fusion1 <- data.frame(sep_Fusion1, check.names = FALSE, stringsAsFactors = FALSE)
clean_Fusion1 <- na.omit(nw.df_Fusion1)
#converts the dataframe to a matrix
m.clean_Fusion1 <- as.matrix(clean_Fusion1)
#converts c(0 or 0) to 0
m.clean_Fusion1[m.clean_Fusion1 =="c(0"] <- 0
m.clean_Fusion1[m.clean_Fusion1 =="c(1"] <- 1
m.clean_Fusion1[m.clean_Fusion1 =="0)"] <- 0
#converts nonzeros to 1
m.clean_Fusion1[m.clean_Fusion1 =="1)"] <- 1
m.clean_Fusion1[m.clean_Fusion1 ==2] <- 1
m.clean_Fusion1[m.clean_Fusion1 =="c(2"] <- 1
m.clean_Fusion1[m.clean_Fusion1 =="2)"] <- 1
m.clean_Fusion1[m.clean_Fusion1 ==3] <- 1 
m.clean_Fusion1[m.clean_Fusion1 =="c(3"] <- 1
m.clean_Fusion1[m.clean_Fusion1 =="3)"] <- 1
m.clean_Fusion1[m.clean_Fusion1 ==4] <- 1 
m.clean_Fusion1[m.clean_Fusion1 =="c(4"] <- 1
m.clean_Fusion1[m.clean_Fusion1 =="4)"] <- 1
m.clean_Fusion1[m.clean_Fusion1 ==5] <- 1 
m.clean_Fusion1[m.clean_Fusion1 =="c(5"] <- 1
m.clean_Fusion1[m.clean_Fusion1 =="5)"] <- 1
m.clean_Fusion1[m.clean_Fusion1 ==6] <- 1 
m.clean_Fusion1[m.clean_Fusion1 =="c(6"] <- 1
m.clean_Fusion1[m.clean_Fusion1 =="6)"] <- 1
m.clean_Fusion1[m.clean_Fusion1 ==7] <- 1 
m.clean_Fusion1[m.clean_Fusion1 =="c(7"] <- 1
m.clean_Fusion1[m.clean_Fusion1 =="7)"] <- 1
m.clean_Fusion1[m.clean_Fusion1 ==8] <- 1
m.clean_Fusion1[m.clean_Fusion1 =="c(8"] <- 1
m.clean_Fusion1[m.clean_Fusion1 =="8)"] <- 1
m.clean_Fusion1[m.clean_Fusion1 ==9] <- 1 
m.clean_Fusion1[m.clean_Fusion1 =="c(9"] <- 1
m.clean_Fusion1[m.clean_Fusion1 =="9)"] <- 1

```

```{r format}
#id.lgf is a temporary variable to store the dataframe of Fusion, Gain, Loss
id.lgf <- data.frame(loss = m.clean_Loss1, gain = m.clean_Gain1, fusion = m.clean_Fusion1, stringsAsFactors = FALSE)
id.lgf$Clones <- clone 
id.lgf$ID <- ID
#id.lgf <- merge(id.lgf,id.index,by="loss.ID",all=TRUE)#merges the id.index to the id.lgf dataframe creates NA since id.index has ids that were previously removed
#id.lgf <- na.omit(id.lgf)#removes any NA data made when merging the id.index
#id.lgf$ID <- paste(id.lgf$loss.ID,id.lgf$index,sep = "_") #produces a new column with the number of replicants 

#id.lgf1<- id.lgf %>% dplyr::select(2751,2750,everything())#rearranges the columns to have the id, and index to the left
id.lgf1<- id.lgf %>% dplyr::select(2750,2749,everything())#rearranges the columns to have the id, and index to the left

#id.lgf1[2:3] <- NULL #removes old id and index

write.csv(id.lgf1, file = "test.csv", row.names = FALSE)
gc()# built in function to free up excess memory
```










