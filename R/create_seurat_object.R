########################################################################################
# scRNA-seq analysis pipeline
# 2 create seurat object
# R 4.1.3; Seurat 4.1.0
# This script reads in raw data from rawData folder
# Single sample mode: barcodes, features, and matrix files are in rawData folder
# Multi-sample mode: bardoces, features, and matrix are organized by samples in rawData
#########################################################################################
library(Seurat)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
inputDir <- args[1]
output <- args[2]
gseName <- args[3]


fs <- list.dirs(inputDir, recursive = F)
samples <- unique(fs)

sceList <- lapply(samples, function(pro) {
    print(paste0("Creating seurat object for sample ", pro))
    CreateSeuratObject(counts = Read10X(pro), project = basename(pro))
})

if (length(sceList) > 1) {
    sce.big <- merge(sceList[[1]],
        y = unlist(sceList)[2:length(samples)],
        add.cell.ids = samples,
        project = gseName
    )
} else {
    sce.big <- sceList[[1]]
}
saveRDS(sce.big, file = output)
