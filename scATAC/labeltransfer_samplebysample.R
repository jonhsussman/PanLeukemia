library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(RColorBrewer)
library(R.utils)
library(ggplot2)
library(readxl)
library(GenomeInfoDb)
library(dplyr)

args = commandArgs(T)
patientName = args[1]

projected.population <- fread("PanLeukemia/scRNA/Metadata/All96Leukemia_proj2healthyCells_to_2021_7_26_full.ref.gini.genes.cc.removed.txt")

projected.population.sample <- projected.population[patient.sample == patientName,]
projected.population.sample <- projected.population.sample[,c("patient_bc", "ProjectedCellType.New", "ProjectedCellType.Binned.New", "CellType")]
rownames(projected.population.sample) <- projected.population.sample$patient_bc

RNA.seurat.path <- "PanLeukemia/scRNA/SeuratObj/BySamples/"
ATAC.mtx.path <- "PanLeukemia/scATAC/scATACpro_Output/"
Plot.path <- "PanLeukemia/scATAC/Figures/Patients/"
seurat.path <- "PanLeukemia/scATAC/SeuratObj/Patients/"
meta.path <- "PanLeukemia/scATAC/Metadata/"
  
seurat.rna <- readRDS(paste0(RNA.seurat.path, patientName, "scRNA.rds"))
seurat.rna <- AddMetaData(seurat.rna, metadata = projected.population.sample)

seurat.atac <- readRDS(paste0(seurat.path, patientName, "_scATAC.rds"))


DefaultAssay(seurat.atac) <- 'GAS'

transfer.anchors <- FindTransferAnchors(
  reference = seurat.rna,
  query = seurat.atac,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = seurat.rna$ProjectedCellType.New,
  weight.reduction = seurat.atac[['lsi']],
  dims = 2:30
)

predicted.labels <- subset(predicted.labels, select = c('predicted.id'))

names(predicted.labels) <- 'labeltransfered.projectedCtype'

seurat.atac <- AddMetaData(object = seurat.atac, metadata = predicted.labels)


saveRDS(seurat.atac, file = paste0(seurat.path, patientName, "_scATAC_withProjectedCellType_labeltransfer.rds"))

meta.final <- seurat.atac@meta.data[,c("nCount_ATAC", "frac_mito", "frac_peak", "predicted.id", "labeltransfered.projectedCtype")]
write.table(meta.final, file = paste0(meta.path, patientName, "_scATAC_withProjectedCellType_labeltransfer.txt"),
            quote = F, sep = "\t")

p3 <- DimPlot(seurat.atac, group.by = "labeltransfered.projectedCtype", label = T) + NoLegend() + labs(title = "scATAC")
p4 <- DimPlot(seurat.rna, group.by = "ProjectedCellType.New", label = T) + NoLegend()+ labs(title = "scRNA")

postscript(paste0(Plot.path, patientName, "_afterlabeltransfer.eps"), height = 6, width = 16)
p3+p4
dev.off()
