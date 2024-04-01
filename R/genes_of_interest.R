########################################################################################
# scRNA-seq analysis pipeline
# 6 check expression of genes of interest
# R 4.1.3; Seurat 4.1.0
########################################################################################

library(Seurat)
library(ggplot2)

# Parameters need to be adjusted for each dataset
args <- commandArgs(trailingOnly = TRUE)
input_data <- args[1]
plot_dir <- args[2]

sc.tcell <- readRDS(input_data)

# options(repr.plot.width=15, repr.plot.height=10)
if ("seurat_clusters" %in% names(sc.tcell@meta.data)) {
    DimPlot(sc.tcell, reduction = "umap", group.by = "seurat_clusters")
    ggsave(file.path(plot_dir, "6_T_cell_umap_seurat_clusters.png"), width = 15, height = 10, units = "in")
}

if ("subject" %in% names(sc.tcell@meta.data)) {
    DimPlot(sc.tcell, reduction = "umap", group.by = "subject")
    ggsave(file.path(plot_dir, "6_T_cell_umap_subject.png"), width = 15, height = 10, units = "in")
}

if ("status" %in% names(sc.tcell@meta.data)) {
    DimPlot(sc.tcell, reduction = "umap", group.by = "status")
    ggsave(file.path(plot_dir, "6_T_cell_umap_status.png"), width = 15, height = 10, units = "in")
}

if ("tissue" %in% names(sc.tcell@meta.data)) {
    DimPlot(sc.tcell, reduction = "umap", group.by = "tissue")
    ggsave(file.path(plot_dir, "6_T_cell_umap_tissue.png"), width = 15, height = 10, units = "in")
}

if ("group" %in% names(sc.tcell@meta.data)) {
    DimPlot(sc.tcell, reduction = "umap", group.by = "group")
    ggsave(file.path(plot_dir, "6_T_cell_umap_group.png"), width = 15, height = 10, units = "in")
}

if ("group" %in% names(sc.tcell@meta.data)) {
    DimPlot(sc.tcell, reduction = "umap", group.by = "group")
    ggsave(file.path(plot_dir, "6_T_cell_umap_group.png"), width = 15, height = 10, units = "in")
}

if ("population" %in% names(sc.tcell@meta.data)) {
    DimPlot(sc.tcell, reduction = "umap", group.by = "population")
    ggsave(file.path(plot_dir, "6_T_cell_umap_population.png"), width = 15, height = 10, units = "in")
}


## One carbon metabolism
FeaturePlot(sc.tcell, c("SHMT1", "SHMT2", "MTHFD1", "MTHFD2", "MTHFR", "TYMS"))
ggsave(file.path(plot_dir, "6_T_cell_umap_one_carbon.png"), width = 10, height = 10, units = "in")

## TCA cycle
FeaturePlot(sc.tcell, c("ACO1", "IDH1", "OGDH", "SUCLA2", "SDHA", "FH", "MDH1", "CS", "ACLY", "LDHA"))
ggsave(file.path(plot_dir, "6_T_cell_umap_TCA.png"), width = 10, height = 10, units = "in")
