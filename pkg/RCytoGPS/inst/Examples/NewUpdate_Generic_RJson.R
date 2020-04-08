# Load the package required to read JSON files.
if (!require("rjson")) {
  install.packages("rjson")
  library("rjson")
}
if (!require("splitstackshape")){
  install.packages("splitstackshape")
  library("splitstackshape")
}
if (!require("dplyr")) {
  install.packages("dplyr")
  library("dplyr")
}
if (!require("jsonlite")) {
  install.packages("jsonlite")
  library("jsonlite")
}
wd <<- choose.dir()
setwd(wd)
message(sprintf("Current working dir: %s\n", wd))


filenames <- list.files(".", pattern="*.json", full.names=TRUE) # this should give you a character vector, with each file name represented by an entry

#doesn't work if you load in more than one json file
myJSON <- lapply(filenames, function(x) rjson::fromJSON(file=filenames)) # a list in which each element is one of your original JSON files


Status <- sapply(myJSON[[1]]$output, "[", "status") #pulls out all the status object list


clone <- unlist(lapply(myJSON[[1]]$output, function(x){ #Function to parse all the loss data from parsing_result
  if((x$status)=="Success"){ 
    1:length(lapply(x$parsing_result, "[[", "loss_gain_fusion_computing"))
  }
}))

rnames <- function(lst, item) {
  f <- function(ll, inds) {
    if ((ii <- match(item, names(ll), FALSE)))
      list(inds=c(inds, ii), len=length(ll[[ii]]))
    else if (all(is.atomic(unlist(ll, FALSE))) || !is.list(ll))
      NULL
    else
      lapply(seq_along(ll), function(i) f(ll[[i]], inds=c(inds, i)))
  }
  unlist(f(lst, NULL))
}


bands <- myJSON[[1]]$iscn2016_bands#pulls bands from top level

lgf <- lapply(myJSON[[1]]$output, function(x){ #Function to parse all the loss data from parsing_result
  if((x$status)=="Success"){ 
    lgf <- lapply(x$parsing_result, "[[","loss_gain_fusion_computing")
  }
})

df.ID <- as.data.frame(matrix(rnames(lgf, "loss"), ncol = 4, byrow = TRUE))
id <- df.ID$V1
loss1 <- sapply(lgf, function(x){
  lapply(x, "[[","loss")
})

gain1 <- sapply(lgf, function(x){
  lapply(x, "[[","gain")
})

fusion1 <- sapply(lgf, function(x){
  lapply(x, "[[","fusion")
})

#seperate the loss, gains, fusion into matrix dataframes and then concatinate them back together
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

#formating the LGF dataframe
df.lgf <- data.frame(loss = m.clean_Loss1, gain = m.clean_Gain1, fusion = m.clean_Fusion1, stringsAsFactors = FALSE)
df.lgf$Clones <- clone
df.lgf$ID <- id

clone.lgf1<- df.lgf %>% dplyr::select(2750,2749,everything())#rearranges the columns to have the id, and index to the left


write.csv(clone.lgf1, file = "NewTest.csv", row.names = FALSE)
gc()# built in function to free up excess memory











