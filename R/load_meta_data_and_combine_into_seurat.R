########################################################################################
# scRNA-seq analysis pipeline
# load meta data from file and combine it with seurat obj
# R 4.1.3;
# by Zhu Liu
# This script load meta data from local file and combine it with the seurat object.
# Recommand to use tibble format for meta data, see https://tibble.tidyverse.org/
#########################################################################################

## TODO: load meta data ##
## TODO: write meta data with seurat object ##
args <- commandArgs(trailingOnly = TRUE)
meta_data <- args[1]
seruatFile <- args[2]
output <- args[3]
gseName <- args[4]


library(Seurat)
library(R.utils)
sc <- readRDS(seruatFile)

df <- read.csv(meta_data, header = TRUE, sep = " ")
for (num in 1:length(df$subject)) {
    if (grepl("Ileal", df$subject[num], fixed = TRUE)) {
        df$subject[num] <- paste0(strsplit(df$subject[num], " ")[[1]][1], " ", strsplit(df$subject[num], " ")[[1]][2])
    }
    if (grepl("PBMC", df$subject[num], fixed = TRUE)) {
        df$subject[num] <- "PBMC"
    }
}

sc@meta.data$subject <- sc@meta.data$orig.ident # lapply(sc@meta.data$orig.ident, function(X) strsplit(X,'_')[[1]][1])
sc@meta.data$subject <- as.factor(unlist(sc@meta.data$subject))
sc@meta.data$status <- lapply(sc@meta.data$subject, function(X) df$condition[which(df$subject == X)])
sc@meta.data$status <- as.factor(unlist(sc@meta.data$status))
saveRDS(sc, file = output)
