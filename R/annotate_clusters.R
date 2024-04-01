########################################################################################
# scRNA-seq analysis pipeline
# 4 annotate clusters
# R 4.1.3; Seurat 4.1.0
# This script read clustered seurat object as the input, and output the automatically
# identified cell identity results.
# It will overwrite the meta data seurat cluster results.
########################################################################################

## load the packages
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)
library(readxl)
library(EnsDb.Hsapiens.v79)
library(SingleR)
library(celldex)

# install EnsDb.Hsapiens.v79, SingleR, celldex if running for the first time
# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("EnsDb.Hsapiens.v79")
# BiocManager::install("SingleR")
# BiocManager::install("celldex")

# Parameters need to be adjusted for each dataset
args <- commandArgs(trailingOnly = TRUE)
input_data <- args[1]
annotated_obj_name <- args[2]
umap_plot <- args[3]

## load the data, the seurat object is named as sc
sc <- readRDS(input_data)

### automatically cell identity identifier
hpca.se <- HumanPrimaryCellAtlasData()
sce_for_SingleR <- GetAssayData(sc, slot = "data")
clusters <- sc@meta.data$seurat_clusters
pred.hesc <- SingleR(
    test = sce_for_SingleR, ref = hpca.se,
    labels = hpca.se$label.fine,
    method = "cluster", clusters = clusters,
    assay.type.test = "logcounts", assay.type.ref = "logcounts"
)


### write the cell identity results to the meta data seurat cluster
seq <- lapply(sc@meta.data$seurat_clusters, function(x) {
    toString(x)
})
for (i in 1:length(rownames(pred.hesc))) {
    seq[which(seq == toString(rownames(pred.hesc)[i]))] <- pred.hesc$labels[i]
}
sc@meta.data$seurat_clusters <- factor(unlist(seq))


### plot the umap
# options(repr.plot.width=15, repr.plot.height=10)
plot1 <- DimPlot(sc, reduction = "umap", group.by = "seurat_clusters")
ggsave(umap_plot, plot = plot1, width = 15, height = 10, units = "in")
saveRDS(sc, file = annotated_obj_name)


# optional:
# install.packages('tidyverse')
#' levels(sc@meta.data$seurat_clusters)
#' library(tidyverse)
#'
#' umap_tx = sc@reductions$umap@cell.embeddings %>%
#' as.data.frame() %>% cbind(tx = sc@meta.data$seurat_clusters)
#'
#' ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2, color=tx)) + geom_point() +
#' scale_color_manual(values=c(#"group1_untreated" = "darkblue",
#'     #'NK_cell'="darkred"
#'     #'T_cell:CD4+'= "darkred",
#'     #'T_cell:CD4+_central_memory'= "darkred",
#'     #'T_cell:CD4+_effector_memory'= "darkred",
#'     #'T_cell:CD4+_Naive'= "darkred"
#'     #'T_cell:CD8+'= "darkred"
#'     #'B_cell:Memory'= "darkred",
#'     'Monocyte'= "darkred",
#'     'Monocyte:CD14+'= "darkred",
#'     'Monocyte:CD16-'= "darkred",
#'     'Monocyte:CD16+'= "darkred",
#'     'Monocyte:leukotriene_D4'= "darkred"
#'     #'T_cell:gamma-delta'= "darkred"
#'                            ))
