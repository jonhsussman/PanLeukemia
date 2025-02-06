source('scDataAnalysis_Utilities.R')
source('dsAnalysis_utilities.R')
library(Seurat)
library(Signac)
library(harmony)
library(ggplot2)
library(parallel)
library(tidyr)
library(dplyr)
library(plyr)
library(readxl)
library(data.table)
library(chromVAR)
library(RColorBrewer)
library(TFBSTools)
library(ggprism)
library(patchwork)
library(magrittr)
library(viridis)
library(DESeq2)


`%notin%` = Negate(`%in%`)

seurat.rna.hspc <- readRDS("HSPC_Integration/seurat_allHSPC_fourSubtypes.rds")
df.cn <- data.table(table(seurat.rna.hspc$patient.sample))
sample.excluded <- df.cn[df.cn$N < 30,]$V1
seurat.rna.hspc <- subset(seurat.rna.hspc, subset = (patient.sample %notin% sample.excluded))

pan.hspc.list <- SplitObject(seurat.rna.hspc, split.by = "patient.sample")

pan.hspc.list <- lapply(X = pan.hspc.list, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = 'vst', nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = pan.hspc.list)


pan.hspc.list <- lapply(X = pan.hspc.list, FUN = function(x){
  x <- ScaleData(x, features = features, verbose = F)
  x <- RunPCA(x, features = features, verbose = F, npcs = 30)
})


pan.hspc.anchors <- FindIntegrationAnchors(object.list = pan.hspc.list, 
                                           anchor.features = features, reduction = "rpca")

pan.hspc.combined <- IntegrateData(anchorset = pan.hspc.anchors, k.weight = 30)


DefaultAssay(pan.hspc.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
pan.hspc.combined <- ScaleData(pan.hspc.combined, verbose = FALSE)
pan.hspc.combined <- RunPCA(pan.hspc.combined, npcs = 30, verbose = FALSE)
pan.hspc.combined <- RunUMAP(pan.hspc.combined, reduction = "pca", dims = 1:30)
pan.hspc.combined <- FindNeighbors(pan.hspc.combined, reduction = "pca", dims = 1:30)
pan.hspc.combined <- FindClusters(pan.hspc.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(pan.hspc.combined, reduction = "umap", group.by = "patient.sample")

p1 <- DimPlot(pan.hspc.combined, reduction = "umap", group.by = "Subtype")


saveRDS(pan.hspc.combined, file = "HSPC_Integration/seurat_rpca_HSPC_fourSubtypes.rds")



