if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("GEOquery")

install.packages("Seurat")
install.packages("stringr")
install.packages("readxl")
install.packages("tidyverse")
BiocManager::install("EnsDb.Hsapiens.v79")
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("pathview")
BiocManager::install("gage")
