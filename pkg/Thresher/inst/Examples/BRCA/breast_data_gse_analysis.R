
############################################################################################
### Consider general cancer profiling data in GEO

## load libraries

library(PCDimension)
library(Thresher)
library(MASS)
library(gdata)
source(file.path("..", "ExistingMethods", "NbClust.txt"))
# library(GEOquery)

# load workspace file that contain all data and analysis results
brcapath <- file.path("scratch", "breast_data_gse_analysis.Rda")
if (file.exists(brcapath)) {
  load(brcapath)
} else { # run the full analysis -- this is long


## read gene list

cid <- read.csv("data/Perou_intrinsic_gene.csv")
perou_annotation <- read.xls("data/Hu-et-al-Intrinsic-List-with-Annotation.xls",
    sheet = 1, header = TRUE)
gene_name <- sapply(1:nrow(cid), function(x) {paste(perou_annotation$Symbol[
    gsub(" ", "", paste(perou_annotation$UGCluster), fixed = TRUE) %in% paste(cid$CLID[x])][1])} )

# gene_name
# cid$CLID[which(gene_name == "")]

gene_ne_name <- gene_name[-which(gene_name == "")]  # 301


############################################################################################
## extract datasets and save as .Rda files

gse_data_list <- c("GSE1992", "GSE2607", "GSE2741", "GSE3143", "GSE3155", "GSE3193", "GSE4611",
    "GSE7422", "GSE10810", "GSE10885", "GSE10886", "GSE12622", "GSE17650", "GSE18539", "GSE19177",
    "GSE19783", "GSE20711", "GSE21921", "GSE22093", "GSE23988", "GSE29431", "GSE37145", "GSE37751",
    "GSE39004", "GSE40115", "GSE43358", "GSE45255", "GSE45827", "GSE46184", "GSE46928", "GSE49481",
    "GSE50939", "GSE53031", "GSE56493", "GSE60785")

## get column names for genes in each datast

gene_ne_name <- gene_name[-which(gene_name == "")]  
load("data/gene_label.Rda")
gene_id_list <- list()
gene_id_list2 <- list()
gene_count_list <- list()
gene_num_list <- list()
gene_name_list <- list()
gene_common_list <- list()
gse_rawdata_list <- list()
gse_data <- list()
for (i in 1:length(gse_data_list)) {
    fi <- paste("RData/GEO_Breast_Data/", gse_data_list[i], ".Rda", sep = '')
    load(fi)
    gene_id_list[[i]] <- list()
    gene_count_list[[i]] <- list()
    gene_num_list[[i]] <- rep(NA, length(gse_i))
    gene_name_list[[i]] <- list()
    gene_id_list2[[i]] <- list()
    gse_rawdata_list[[i]] <- list()
    gse_data[[i]] <- NULL
    for (k in 1:length(gse_i)) {
      if (!is.na(gene_label_list[[i]][k])) {
        gene_id_list[[i]][[k]] <- lapply(1:length(gene_ne_name), function(x) {which(fData(gse_i[[k]])[, 
          gene_label_list[[i]][k]] %in% gene_ne_name[x])} )
        gene_count_list[[i]][[k]] <- sapply(gene_id_list[[i]][[k]], length)
        gene_num_list[[i]][k] <- sum(gene_count_list[[i]][[k]] > 0)
        gene_name_list[[i]][[k]] <- gene_ne_name[gene_count_list[[i]][[k]] > 0]
      } else {
        gene_id_list[[i]][[k]] <- NA
        gene_count_list[[i]][[k]] <- NA
        gene_num_list[[i]][k] <- NA
        gene_name_list[[i]][[k]] <- NA
      }
    }
    gene_common_list[[i]] <- gene_ne_name
    for (k in 1:length(gse_i)) {
      if (!is.na(gene_label_list[[i]][k])) {
        gene_common_list[[i]] <- intersect(gene_common_list[[i]], gene_name_list[[i]][[k]])
      }
    }
    gse_data_i <- NULL
    for (k in 1:length(gse_i)) {
      if (!is.na(gene_label_list[[i]][k])) {
        gene_id_list2[[i]][[k]] <- sapply(1:length(gene_common_list[[i]]), function(x) {
          which(fData(gse_i[[k]])[, gene_label_list[[i]][k]] %in% gene_common_list[[i]][x])[1]} )
        gse_rawdata_list[[i]][[k]] <- exprs(gse_i[[k]])[gene_id_list2[[i]][[k]], ]
        rownames(gse_rawdata_list[[i]][[k]]) <- gene_common_list[[i]]
        # different GPL may have different normalization 
        gse_data_i <- cbind(gse_data_i, t(scale(t(gse_rawdata_list[[i]][[k]])))) 
      }
    }
    gse_data[[i]] <- gse_data_i 

    rm(gse_i)
}

# lapply(gse_data, dim)

## only look at the data set of appropriate size

ind <- c(1:7, 9:12, 15:22, 24:30, 32:35)


############################################################################################
### perform Thresher anlysis for each dataset

# list criteria in PCDimension to compute the number of PCs
f <- makeAgCpmFun("Exponential")
agfuns <- list(twice = agDimTwiceMean, specc = agDimSpectral, km = agDimKmeans, 
               km3 = agDimKmeans3, tt = agDimTtest, tt2 = agDimTtest2,
               cpt = agDimCPT, cpm = f)

# list all indices in NbClust package
indfun <- c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", 
           "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", 
           "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", 
           "gamma", "gplus", "tau", "dunn", "sdindex", "sdbw")
minnc <- 2
maxnc <- 10

pcdim <- rep(0, length(ind))
noutlier <- rep(0, length(ind))
nsample <- rep(0, length(ind))

thresher.nc <- rep(0, length(ind))
nb.nc.1ist.1 <- list()
outlierlist.1 <- list()
sdir <- "Plots_and_Results/breast_data_gse"

## use TwiceMean criterion to get the number of clusters for each dataset

for (k in 1:length(ind)) {
  
  j <- ind[k]
  sdir_j <- paste(sdir, "/", gse_data_list[j], sep = '')
  if(!file.exists(sdir_j)) dir.create(sdir_j)

  sink(paste(sdir_j, "/", gse_data_list[j], ".txt", sep=''))

  # clean data since it may have missing values
  data <- gse_data[[j]]
  nsample[k] <- ncol(data)
  data[is.na(data)] <- 0

  # use Thresher to estimate number of clusters

  spca1 <- SamplePCA(t(data))
  # spca1 <- SamplePCA(scale(t(data)))
  obj1 <- AuerGervini(spca1)
  
  print(compareAgDimMethods(obj1, agfuns))
 
  png(paste(sdir_j, "/", gse_data_list[j], "_ag.png", sep=''))
  plot(obj1, agfuns, ylim=c(0, 30))
  dev.off()

  thresh1 <- Thresher(as.matrix(data), method = "auer.gervini", scale = FALSE, agfun = 
    agDimTwiceMean)
  # thresh1 <- Thresher(t(scale(t(data))), method = "auer.gervini", scale = FALSE, agfun = 
  #   agDimTwiceMean)
  noutlier[k] <- sum(thresh1@delta <= 0.3)
  outlierlist.1[[k]] <- colnames(data)[thresh1@delta <= 0.3]
  cat("PC Dimension is", thresh1@pcdim, "\n") 
  cat("Number of outliers is", sum(thresh1@delta <= 0.3), "\n") 

  reaper1 <- Reaper(thresh1, cutoff = 0.3, useLoadings = TRUE, method = "auer.gervini")
  cat("Number of Groups is", reaper1@nGroups, "\n")
  print(reaper1@bic)

  pcdim[k] <- thresh1@pcdim
  thresher.nc[k] <- reaper1@nGroups

  png(paste(sdir_j, "/", gse_data_list[j], "_thresh.png", sep=''))
  plot(thresh1)
  dev.off()
  png(paste(sdir_j, "/", gse_data_list[j], "_reaper.png", sep=''))
  plot(reaper1)
  dev.off()

  ## use NbClust package indices to perform analysis 

  # using NbClust indices to estimate number of clusters
  nb.nc.1ist.1[[k]] <- rep(0, length(indfun))
  for (i in 1:length(indfun)) {
    fit <- try(NbClust(t(data), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
          method = "ward.D2", index = indfun[i]))
    if (class(fit) != "try-error") { 
      nb.nc.1ist.1[[k]][i] <- as.numeric(fit$Best.nc[1])
      cat("For index ", indfun[i], ", the number of clusters is ", as.numeric(fit$Best.nc[1]), 
        sep = '')
      cat("\n")
    } else {
      nb.nc.1ist.1[[k]][i] <- NA
      cat("No results return for index \n") 
    }
  }

  tabname <- paste(sdir_j, "/", gse_data_list[j], "_nb.nc.csv", sep='')
  write.table(cbind(index = indfun, number = nb.nc.1ist.1[[k]]), tabname, sep = ",", 
    col.names = NA, row.names = TRUE)

  sink()

}

# obtain the results of clustering for all datasets
thresher.mat.1 <- data.frame(cbind(GSE = gse_data_list[ind], PCDimension = pcdim, Thresher_NC = thresher.nc))
outlier_mat <- cbind(sample_size = nsample, outlier_number = noutlier, outlier_percentage = 100*noutlier/nsample)
rownames(outlier_mat) <- gse_data_list[ind]
outlier_mat_v <- outlier_mat[gse_data_list2, ]

# write.table(thresher.mat.1, "GSE_TwiceMean_result.csv", sep = ",", col.names = NA, row.names = TRUE)
# write.table(outlier_mat_v, "GSE_TwiceMean_Outlier_Summary.csv", sep = ",", col.names = NA, row.names = TRUE)


## use CPT criterion to get the number of clusters for each dataset

pcdim2 <- rep(0, length(ind))
noutlier2 <- rep(0, length(ind))
thresher.nc2 <- rep(0, length(ind))
nb.nc.1ist.2 <- list()
outlierlist.2 <- list()

for (k in 1:length(ind)) {
  
  j <- ind[k]
  sdir_j <- paste(sdir, "/", gse_data_list[j], sep = '')
  if(!file.exists(sdir_j)) dir.create(sdir_j)

  sink(paste(sdir_j, "/", gse_data_list[j], ".txt", sep=''))

  # clean data since it may have missing values
  data <- gse_data[[j]]
  data[is.na(data)] <- 0

  # use Thresher to estimate number of clusters

  spca1 <- SamplePCA(t(data))
  # spca1 <- SamplePCA(scale(t(data)))
  obj1 <- AuerGervini(spca1)

  print(compareAgDimMethods(obj1, agfuns))
 
  png(paste(sdir_j, "/", gse_data_list[j], "_ag.png", sep=''))
  plot(obj1, agfuns, ylim=c(0, 30))
  dev.off()

  thresh1 <- Thresher(as.matrix(data), method = "auer.gervini", scale = FALSE, agfun = 
    agDimCPT)
  # thresh1 <- Thresher(t(scale(t(data))), method = "auer.gervini", scale = FALSE, agfun = 
  #   agDimCPT)
  noutlier2[k] <- sum(thresh1@delta <= 0.3)
  outlierlist.2[[k]] <- colnames(data)[thresh1@delta <= 0.3]
  cat("PC Dimension is", thresh1@pcdim, "\n") 
  cat("Number of outliers is", sum(thresh1@delta <= 0.3), "\n") 

  reaper1 <- Reaper(thresh1, cutoff = 0.3, useLoadings = TRUE, method = "auer.gervini")
  cat("Number of Groups is", reaper1@nGroups, "\n")
  print(reaper1@bic)

  pcdim2[k] <- thresh1@pcdim
  thresher.nc2[k] <- reaper1@nGroups

  png(paste(sdir_j, "/", gse_data_list[j], "_thresh.png", sep=''))
  plot(thresh1)
  dev.off()
  png(paste(sdir_j, "/", gse_data_list[j], "_reaper.png", sep=''))
  plot(reaper1)
  dev.off()

  ## use NbClust package indices to perform analysis 

  # using NbClust indices to estimate number of clusters
  nb.nc.1ist.2[[k]] <- rep(0, length(indfun))
  for (i in 1:length(indfun)) {
    fit <- try(NbClust(t(data), distance = "euclidean", min.nc = minnc, max.nc = maxnc,
          method = "ward.D2", index = indfun[i]))
    if (class(fit) != "try-error") { 
      nb.nc.1ist.2[[k]][i] <- as.numeric(fit$Best.nc[1])
      cat("For index ", indfun[i], ", the number of clusters is ", as.numeric(fit$Best.nc[1]), 
        sep = '')
      cat("\n")
    } else {
      nb.nc.1ist.2[[k]][i] <- NA
      cat("No results return for index \n")
    }
  }

  tabname <- paste(sdir_j, "/", gse_data_list[j], "_nb.nc.csv", sep='')
  write.table(cbind(index = indfun, number = nb.nc.1ist.2[[k]]), tabname, sep = ",", 
    col.names = NA, row.names = TRUE)

  sink()

}

thresher.mat.2 <- data.frame(cbind(GSE = gse_data_list[ind], PCDimension = pcdim2, Thresher_NC = thresher.nc2))

# write.table(thresher.mat.2, "GSE_CPT_result.csv", sep = ",", col.names = NA, row.names = TRUE)


## save the result of NbClust indices in a matrix

# use matrix to save number of clusters based on NbClust indices (nb.nc.1ist.1 and nb.nc.1ist.2 are 
# the same)
nb.nc.mat.1 <- matrix(unlist(nb.nc.1ist.1), nrow = length(ind), ncol = length(indfun), byrow = TRUE)
# nb.nc.mat.2 <- matrix(unlist(nb.nc.1ist.2), nrow = length(ind), ncol = length(indfun), byrow = TRUE)
rownames(nb.nc.mat.1) <- gse_data_list[ind]
colnames(nb.nc.mat.1) <- indfun

# write.table(nb.nc.mat.1, "GSE_NbClust_result.csv", sep = ",", col.names = NA, row.names = TRUE)
# write.table(nb.nc.mat.2, "GSE_NbClust_result_2.csv", sep = ",", col.names = NA, row.names = TRUE)

top10 <- c("trcovw", "tracew", "ratkowsky", "mcclain", "ptbiserial", "tau", "sdindex", "kl", "ccc",
    "hartigan")

# use matrix to save number of clusters for the top 10 NbClust indices and Thresher with 2 best criteria
res <- cbind(nb.nc.mat.1[, top10], Thresher_TwiceMean = thresher.nc, Thresher_CPT = thresher.nc2)
# write.table(res, "GSE_Comparison_result.csv", sep = ",", col.names = NA, row.names = TRUE)

  save.image(file=f)
} # end the full analysis
rm(brcapath)
ls()


# plot the number of clusters
par(mfrow = c(4, 3))
for (i in 1:12) {
    hist(res[, i], xlim = c(0.5, 10.5), xlab = "Number of clusters", breaks = seq(0.5, 10.5, by = 1), 
      main = paste("Histogram of clusters #:", colnames(res)[i]))
    # abline(v = 6, lty = 2, col = 3)
}

par(mfrow = c(3, 5))
for (i in 1:15) {
    hist(nb.nc.mat.1[, i], xlim = c(0.5, 10.5), xlab = "Number of clusters", breaks = seq(0.5, 10.5, 
      by = 1), main = paste("Histogram of clusters #:", colnames(nb.nc.mat.1)[i]))
    # abline(v = 6, lty = 2, col = 3)
}

par(mfrow = c(3, 5))
for (i in 15+1:13) {
    hist(nb.nc.mat.1[, i], xlim = c(0.5, 10.5), xlab = "Number of clusters", breaks = seq(0.5, 10.5, 
      by = 1), main = paste("Histogram of clusters #:", colnames(nb.nc.mat.1)[i]))
    # abline(v = 6, lty = 2, col = 3)
}

# use matrix to save number of clusters for the top 11 NbClust indices and Thresher with TwiceMean criteria
top11 <- c("trcovw", "tracew", "ratkowsky", "mcclain", "ptbiserial", "tau", "sdindex", "kl", "ccc",
    "hartigan", "sdbw")
res2 <- cbind(nb.nc.mat.1[, top11], Thresher_TwiceMean = thresher.nc)

png("hist.png", width = 2*960, height = 2*1280, units = "px", pointsize = 32)
# jpeg("hist.jpeg", width = 2*960, height = 2*1280, units = "px", pointsize = 32, quality = 100)
par(mfrow = c(4, 3))
for (i in 1:12) {
    hist(res2[, i], xlim = c(0.5, 10.5), xlab = "Number of clusters", breaks = seq(0.5, 10.5, by = 1), 
      main = paste("Histogram for method :", colnames(res2)[i]))
    # abline(v = 6, lty = 2, col = 3)
}
dev.off()

# remove 5 datasets that has small samples sizes (<= about 60)
gse_data_list2 <- setdiff(gse_data_list[ind], c("GSE3155", "GSE3193", "GSE10886", "GSE23988", "GSE46928"))


# plot the results (histograms)

png("hist2.png", width = 2*960, height = 2*1280, units = "px", pointsize = 32)
# jpeg("hist2.jpeg", width = 2*960, height = 2*1280, units = "px", pointsize = 32, quality = 100)
par(mfrow = c(4, 3))
for (i in 1:12) {
    hist(res2[gse_data_list2, i], xlim = c(0.5, 10.5), xlab = "Number of clusters", breaks = seq(0.5, 
      10.5, by = 1), main = paste("Histogram for method :", colnames(res2)[i]))
    # abline(v = 6, lty = 2, col = 3)
}
dev.off()

pdf("hist2.pdf")
par(mfrow = c(4, 3))
par(mar = c(2, 2, 2, 2))
for (i in 1:12) {
    hist(res2[gse_data_list2, i], xlim = c(0.5, 10.5), xlab = "Number of clusters", breaks = seq(0.5, 
      10.5, by = 1), main = paste("Histogram for", colnames(res2)[i]))
    # abline(v = 6, lty = 2, col = 3)
}
dev.off()




