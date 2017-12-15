
### Download GSE data by using GEOquery package

library(GEOquery)

gse_data_list <- c("GSE1992", "GSE2607", "GSE2741", "GSE3143", "GSE3155", "GSE3193", "GSE4611",
    "GSE7422", "GSE10810", "GSE10885", "GSE10886", "GSE12622", "GSE17650", "GSE18539", "GSE19177",
    "GSE19783", "GSE20711", "GSE21921", "GSE22093", "GSE23988", "GSE29431", "GSE37145", "GSE37751",
    "GSE39004", "GSE40115", "GSE43358", "GSE45255", "GSE45827", "GSE46184", "GSE46928", "GSE49481",
    "GSE50939", "GSE53031", "GSE56493", "GSE60785")

dir1 <- "RData"
if (!file.exists(dir1)) dir.create(dir1)
dir2 <- "RData//GEO_Breast_Data"
if (!file.exists(dir2)) dir.create(dir2)

for (i in 1:length(gse_data_list)) {
    gse_num_i <- gse_data_list[i]
    gse_i <- getGEO(gse_num_i)
    save(gse_i, file = paste("RData/GEO_Breast_Data/", gse_num_i, ".Rda", sep = ''))
}

