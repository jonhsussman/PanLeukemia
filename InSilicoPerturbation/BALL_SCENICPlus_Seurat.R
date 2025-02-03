library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(anndata)
library(Signac)
library(Matrix)
library(ggplot2)
library(dplyr)
library(SeuratDisk)

setwd("/mnt/isilon/tan_lab/sussmanj/Temp/PanLeukemia/SCENICPlus")

BALL_metacells <- readRDS("/mnt/isilon/tan_lab/sussmanj/Temp/PanLeukemia/SCENICPlus/BALL_Metacells_Projected.rds")

#SCENIC+ Results
path = "/mnt/isilon/tan_lab/sussmanj/Temp/PanLeukemia/SCENICPlus/SCENICPlus_Pipeline/B_ALL/Output/"
Convert(paste0(path, "eGRN_AUC.h5ad"), dest = "h5seurat", overwrite = TRUE)
scenicplus <- LoadH5Seurat(paste0(path, "eGRN_AUC.h5seurat"), meta.data = FALSE, misc = FALSE)
rownames(scenicplus)
suffix <- "___cisTopic"
cellnames_correct = colnames(BALL_metacells)

eGRN = GetAssayData(scenicplus, assay = "RNA", slot = "data")
colnames(eGRN) = sub("_+cisTopic$", "", colnames(eGRN))
eGRN_reordered <- eGRN[, match(cellnames_correct, colnames(eGRN))]
head(colnames(eGRN_reordered))
rownames(eGRN_reordered) <- gsub("(?<!/)-(?=[^/])", "_", rownames(eGRN_reordered), perl = TRUE)
identical(cellnames_correct, colnames(eGRN_reordered))

eGRN_metadata = read.table(paste0(path, "eRegulon_filtered.tsv"), sep = '\t', header = T)
head(eGRN_metadata)
filtered = c(eGRN_metadata$Gene_signature_name, eGRN_metadata$Region_signature_name)
eGRN_filtered = eGRN_reordered[rownames(eGRN_reordered) %in% filtered, ]
dim(eGRN_filtered)
rownames(eGRN_filtered)

BALL_metacells[["eGRN"]] <- CreateAssayObject(data = eGRN_filtered)
DefaultAssay(BALL_metacells) =  "eGRN"
#HOXA9_direct_+/+_(159g)
FeaturePlot(BALL_metacells, features = "HOXA9-direct-+/+-(159g)", reduction = "ref.umap", max.cutoff = 'q99', min.cutoff = 'q1', order = T) + coord_fixed()
VlnPlot(BALL_metacells, features = "HOXA9-direct-+/+-(159g)", group.by = "predicted.celltype_refmap", pt.size = 0)
VlnPlot(BALL_metacells, features = "HOXA9-direct-+/+-(159g)", group.by = "predicted.celltypegeneral_refmap", pt.size = 0)

Idents(BALL_metacells) = "predicted.celltype_refmap"
markers = FindAllMarkers(BALL_metacells, assay = "eGRN", logfc.threshold = 0.001)

saveRDS(BALL_metacells, "BALL_Metacells_Projected_SCENICplus.rds")

BALL_metacells = readRDS("BALL_Metacells_Projected_SCENICplus.rds")
write.csv(BALL_metacells@reductions$ref.umap@cell.embeddings, file='BALL_metacells_refUMAP.csv', quote=F, row.names=T)
DimPlot(BALL_metacells, reduction = "ref.umap", group.by = "predicted.trajectory") + coord_fixed()


