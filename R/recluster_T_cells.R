########################################################################################
# scRNA-seq analysis pipeline
# 5 select T cells and recluster
# R 4.1.3; Seurat 4.1.0
# This script take cell type identified and clustered Seurat object. Then pick out all T cell clusters.
########################################################################################

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)
library(readxl)
library(EnsDb.Hsapiens.v79)
library(SingleR)
library(celldex)

# Parameters need to be adjusted for each dataset
memory.limit(9999999999)
args <- commandArgs(trailingOnly = TRUE)
input_data <- args[1]
t_recluster_obj_name <- args[2]
cell_type_cluster_umap <- args[3]
cell_type <- args[4]

print(cell_type)
sc <- readRDS(input_data)
## find the T cell clusters
subset_cluster_number <- which(lapply(levels(sc@meta.data$seurat_clusters), function(x) {
    grepl(cell_type, x, fixed = TRUE)
}) == TRUE)

### subset only the T cells
sc.tcell <- subset(sc, cells = sc@assays$RNA@data@Dimnames[[2]][
    which(sapply(as.numeric(sc@meta.data$seurat_clusters), function(x) x %in% subset_cluster_number))
])

### get cluster unique

old_cluster_name <- unique(levels(sc.tcell@meta.data$seurat_clusters))
sc.tcell <- RunPCA(sc.tcell, features = VariableFeatures(object = sc.tcell))
sc.tcell <- RunUMAP(object = sc.tcell, reduction = "harmony", dims = 1:20, verbose = F)
sc.tcell <- FindNeighbors(object = sc.tcell, reduction = "harmony", dims = 1:20, verbose = F)
sc.tcell <- FindClusters(object = sc.tcell, resolution = 0.5, verbose = F)
hpca.se <- HumanPrimaryCellAtlasData()
sce_for_SingleR <- GetAssayData(sc.tcell, slot = "data")
clusters <- sc.tcell@meta.data$seurat_clusters
pred.hesc <- SingleR(
    test = sce_for_SingleR, ref = hpca.se,
    labels = hpca.se$label.fine,
    method = "cluster", clusters = clusters,
    assay.type.test = "logcounts", assay.type.ref = "logcounts"
)

### write the cell identity results to the meta data seurat cluster
seq <- lapply(sc.tcell@meta.data$seurat_clusters, function(x) {
    toString(x)
})
for (i in 1:length(rownames(pred.hesc))) {
    seq[which(seq == toString(rownames(pred.hesc)[i]))] <- pred.hesc$labels[i]
}
sc.tcell@meta.data$seurat_clusters <- factor(unlist(seq))

options(repr.plot.width = 15, repr.plot.height = 10)
plot1 <- DimPlot(sc.tcell, reduction = "umap", group.by = "seurat_clusters")

ggsave(cell_type_cluster_umap, plot = plot1, width = 15, height = 10, units = "in")
saveRDS(sc.tcell, file = t_recluster_obj_name)
