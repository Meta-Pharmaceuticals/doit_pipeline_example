########################################################################################
# scRNA-seq analysis pipeline
# 3 run seurat
# R 4.1.3; Seurat 4.1.0
# by Zikun Wang
# This script loads the .rds file saved in Step 2 and run seurat pipeline on the data
# Processed seurat object will be saved
#########################################################################################

# Parameters need to be adjusted for each dataset
args <- commandArgs(trailingOnly = TRUE)
input_data <- args[1]
processed_obj_name <- args[2]
QC_metrics_plot <- args[3]
umap_plot <- args[4]

min_genes <- strtoi(args[5])
max_mt_pct <- strtoi(args[6])

library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)

memory.limit(9999999999)
# load r object from Step 2
sc <- readRDS(input_data)

# Calculate percentage of mitochondria genes
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
sc <- SetIdent(object = sc, value = sc@project.name)
# Visualize QC metrics as a violin plot
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(QC_metrics_plot, plot = plot1 + plot2)


#### parameters in the subset function need to be adjusted for each dataset ####
# Filter out cells with less than 1000 nFeature and more than 25% percent.mt
sc <- subset(sc, subset = nFeature_RNA > min_genes & percent.mt < max_mt_pct)

# normalize and scale the data
sc <- NormalizeData(sc, verbose = F)
sc <- FindVariableFeatures(sc, verbose = T, nfeatures = 3000)
sc <- ScaleData(sc)
# sc <- SCTransform(sc, vars.to.regress = c("percent.mt"),ncells=5000)#conserve.memory = TRUE,

# run Seurat clustering
sc <- RunPCA(sc, features = VariableFeatures(object = sc)) #
# Integration
sc <- RunHarmony(sc,
    group.by.vars = c("subject"),
    reduction = "pca", reduction.save = "harmony"
) # ,  assay.use = "SCT"

sc <- RunUMAP(object = sc, reduction = "harmony", dims = 1:20, verbose = F) # , assay = "SCT"
sc <- FindNeighbors(object = sc, dims = 1:20, verbose = F, reduction = "harmony")
sc <- FindClusters(object = sc, resolution = 0.5, verbose = F)
plot3 <- DimPlot(sc, reduction = "umap")
ggsave(umap_plot, plot = plot3)
saveRDS(sc, file = processed_obj_name)
