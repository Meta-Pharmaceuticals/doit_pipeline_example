########################################################################################
# scRNA-seq analysis pipeline
# download meta data
# R 4.1.3;
# by Zhu Liu
# This script downloads meta data from GSE repository and annotates the seurat object.
# Recommand to use tibble format for meta data, see https://tibble.tidyverse.org/
#########################################################################################

## TODO: convert to cli args ##
args <- commandArgs(trailingOnly = TRUE)
output <- args[1]
gseName <- args[2]


# gse_name <- "GSE116222"
metadata_file_name <- "meta.csv"

library(GEOquery)
library(stringr)
library(tidyverse)

# download metadata from GEO
gse.info <- getGEO(gseName, GSEMatrix = TRUE)

# some GSE returns a list of multiple entries, unpack the list
gse.meta <- NULL
for (info in gse.info) {
    # the meta data is stored in phenoData for each element in the list
    info.meta <- info@phenoData@data
    gse.meta <- rbind(gse.meta, info.meta)
}

metadata <- tibble::rownames_to_column(gse.meta, var = "sample") %>%
    as_tibble() %>%
    dplyr::select("geo_accession", "title", "source_name_ch1") %>%
    dplyr::rename(sample = "geo_accession", subject = "title", condition = "source_name_ch1")

write.table(metadata, output)
