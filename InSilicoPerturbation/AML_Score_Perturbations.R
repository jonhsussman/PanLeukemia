library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(anndata)
library(Signac)
library(Matrix)
library(ggplot2)
library(dplyr)
library(SeuratDisk)
library(dittoSeq)
library(AUCell)
library(ggplot2)
library(gridExtra)

setwd("/mnt/isilon/tan_lab/sussmanj/Temp/PanLeukemia/SCENICPlus")

AML_metacells <- readRDS("/mnt/isilon/tan_lab/sussmanj/Temp/PanLeukemia/SCENICPlus/AML_Metacells_Projected_SCENICplus.rds")
DefaultAssay(AML_metacells) = "RNA"
AML_metacells = JoinLayers(AML_metacells)
normal = as.matrix(GetAssayData(AML_metacells, assay = "RNA", layer = "data"))

dittoBarPlot(AML_metacells, var = "predicted.celltype_refmap", group.by = "orig.ident")
dittoBarPlot(AML_metacells, var = "predicted.celltypegeneral_refmap", group.by = "orig.ident")

healthy_ref = readRDS('/mnt/isilon/tan_lab/sussmanj/Temp/T_ALL/Healthy_Reference_Fixed.rds')

hspc_trn = read.table("/mnt/isilon/tan_lab/chenc6/PanLeukemia/scATAC/Scripts/TRN_Updated_withAdditionalBALL/ExampleTRN_Visualization/Shared_TRN_edges_coreSet.txt", header = T)
hspc_genes = unique(c(hspc_trn$TF, hspc_trn$gene_name))

########################################################3
#####JUNB
########################################################
perturbation_data <- read.table("SCENICPlus_Pipeline/AML/Perturbation/AML_JUNB_perturbation_iter_11_processed.tsv", header = TRUE, row.names = 1, sep = "\t")
perturbation_data = t(perturbation_data)
aligned_perturbation <- perturbation_data[, colnames(AML_metacells)]
perturbation_matrix <- as.matrix(aligned_perturbation)
perturbation_matrix <- pmax(perturbation_matrix, 0)
AML_metacells[["JUNB_Perturbation"]] <- CreateAssay5Object(data = perturbation_matrix)
print(AML_metacells)

perturations = as.matrix(GetAssayData(AML_metacells, assay = "JUNB_Perturbation", layer = "data"))
metadata_perturb = data.frame(Perturbed = c(rep("NoPerturb", dim(AML_metacells)[2]), rep("Perturb", dim(AML_metacells)[2])))
colnames(perturations) = paste0(colnames(perturations), "_perturb")
combined_perturbations = cbind(normal, perturations)
dim(combined_perturbations)
rownames(metadata_perturb) = colnames(combined_perturbations)
combined_perturbations_seurat = CreateSeuratObject(counts = combined_perturbations, assay = "RNA")
combined_perturbations_seurat[["RNA"]] = CreateAssay5Object(data = combined_perturbations)
combined_perturbations_seurat = AddMetaData(combined_perturbations_seurat, metadata = metadata_perturb)
table(combined_perturbations_seurat$Perturbed)
Idents(combined_perturbations_seurat) = "Perturbed"
perturbed_degs = FindMarkers(combined_perturbations_seurat, only.pos = F, ident.1 = "Perturb",
                             min.pct = 0, min.diff.pct = 0, logfc.threshold = 0)
saveRDS(perturbed_degs, "DEGs/AML_JUNB_Perturbation.rds")

########################################################3
#####HOXA9
########################################################
perturbation_data <- read.table("SCENICPlus_Pipeline/AML/Perturbation/AML_HOXA9_perturbation_iter_11_processed.tsv", header = TRUE, row.names = 1, sep = "\t")
perturbation_data = t(perturbation_data)
aligned_perturbation <- perturbation_data[, colnames(AML_metacells)]
perturbation_matrix <- as.matrix(aligned_perturbation)
perturbation_matrix <- pmax(perturbation_matrix, 0)
AML_metacells[["HOXA9_Perturbation"]] <- CreateAssay5Object(data = perturbation_matrix)
print(AML_metacells)

perturations = as.matrix(GetAssayData(AML_metacells, assay = "HOXA9_Perturbation", layer = "data"))
metadata_perturb = data.frame(Perturbed = c(rep("NoPerturb", dim(AML_metacells)[2]), rep("Perturb", dim(AML_metacells)[2])))
colnames(perturations) = paste0(colnames(perturations), "_perturb")
combined_perturbations = cbind(normal, perturations)
dim(combined_perturbations)
rownames(metadata_perturb) = colnames(combined_perturbations)
combined_perturbations_seurat = CreateSeuratObject(counts = combined_perturbations, assay = "RNA")
combined_perturbations_seurat[["RNA"]] = CreateAssay5Object(data = combined_perturbations)
combined_perturbations_seurat = AddMetaData(combined_perturbations_seurat, metadata = metadata_perturb)
table(combined_perturbations_seurat$Perturbed)
Idents(combined_perturbations_seurat) = "Perturbed"
perturbed_degs = FindMarkers(combined_perturbations_seurat, only.pos = F, ident.1 = "Perturb",
                             min.pct = 0, min.diff.pct = 0, logfc.threshold = 0)
saveRDS(perturbed_degs, "DEGs/AML_HOXA9_Perturbation.rds")


########################################################3
#####HOXA5
########################################################
perturbation_data <- read.table("SCENICPlus_Pipeline/AML/Perturbation/AML_HOXA5_perturbation_iter_11_processed.tsv", header = TRUE, row.names = 1, sep = "\t")
perturbation_data = t(perturbation_data)
aligned_perturbation <- perturbation_data[, colnames(AML_metacells)]
perturbation_matrix <- as.matrix(aligned_perturbation)
perturbation_matrix <- pmax(perturbation_matrix, 0)
AML_metacells[["HOXA5_Perturbation"]] <- CreateAssay5Object(data = perturbation_matrix)
print(AML_metacells)

perturations = as.matrix(GetAssayData(AML_metacells, assay = "HOXA5_Perturbation", layer = "data"))
metadata_perturb = data.frame(Perturbed = c(rep("NoPerturb", dim(AML_metacells)[2]), rep("Perturb", dim(AML_metacells)[2])))
colnames(perturations) = paste0(colnames(perturations), "_perturb")
combined_perturbations = cbind(normal, perturations)
dim(combined_perturbations)
rownames(metadata_perturb) = colnames(combined_perturbations)
combined_perturbations_seurat = CreateSeuratObject(counts = combined_perturbations, assay = "RNA")
combined_perturbations_seurat[["RNA"]] = CreateAssay5Object(data = combined_perturbations)
combined_perturbations_seurat = AddMetaData(combined_perturbations_seurat, metadata = metadata_perturb)
table(combined_perturbations_seurat$Perturbed)
Idents(combined_perturbations_seurat) = "Perturbed"
perturbed_degs = FindMarkers(combined_perturbations_seurat, only.pos = F, ident.1 = "Perturb",
                             min.pct = 0, min.diff.pct = 0, logfc.threshold = 0)
saveRDS(perturbed_degs, "DEGs/AML_HOXA5_Perturbation.rds")


########################################################3
#####HOXA3
########################################################
perturbation_data <- read.table("SCENICPlus_Pipeline/AML/Perturbation/AML_HOXA3_perturbation_iter_11_processed.tsv", header = TRUE, row.names = 1, sep = "\t")
perturbation_data = t(perturbation_data)
aligned_perturbation <- perturbation_data[, colnames(AML_metacells)]
perturbation_matrix <- as.matrix(aligned_perturbation)
perturbation_matrix <- pmax(perturbation_matrix, 0)
AML_metacells[["HOXA3_Perturbation"]] <- CreateAssay5Object(data = perturbation_matrix)
print(AML_metacells)

perturations = as.matrix(GetAssayData(AML_metacells, assay = "HOXA3_Perturbation", layer = "data"))
metadata_perturb = data.frame(Perturbed = c(rep("NoPerturb", dim(AML_metacells)[2]), rep("Perturb", dim(AML_metacells)[2])))
colnames(perturations) = paste0(colnames(perturations), "_perturb")
combined_perturbations = cbind(normal, perturations)
dim(combined_perturbations)
rownames(metadata_perturb) = colnames(combined_perturbations)
combined_perturbations_seurat = CreateSeuratObject(counts = combined_perturbations, assay = "RNA")
combined_perturbations_seurat[["RNA"]] = CreateAssay5Object(data = combined_perturbations)
combined_perturbations_seurat = AddMetaData(combined_perturbations_seurat, metadata = metadata_perturb)
table(combined_perturbations_seurat$Perturbed)
Idents(combined_perturbations_seurat) = "Perturbed"
perturbed_degs = FindMarkers(combined_perturbations_seurat, only.pos = F, ident.1 = "Perturb",
                             min.pct = 0, min.diff.pct = 0, logfc.threshold = 0)
saveRDS(perturbed_degs, "DEGs/AML_HOXA3_Perturbation.rds")


########################################################3
#####FOSB
########################################################
perturbation_data <- read.table("SCENICPlus_Pipeline/AML/Perturbation/AML_FOSB_perturbation_iter_11_processed.tsv", header = TRUE, row.names = 1, sep = "\t")
perturbation_data = t(perturbation_data)
aligned_perturbation <- perturbation_data[, colnames(AML_metacells)]
perturbation_matrix <- as.matrix(aligned_perturbation)
perturbation_matrix <- pmax(perturbation_matrix, 0)
AML_metacells[["FOSB_Perturbation"]] <- CreateAssay5Object(data = perturbation_matrix)
print(AML_metacells)

perturations = as.matrix(GetAssayData(AML_metacells, assay = "FOSB_Perturbation", layer = "data"))
metadata_perturb = data.frame(Perturbed = c(rep("NoPerturb", dim(AML_metacells)[2]), rep("Perturb", dim(AML_metacells)[2])))
colnames(perturations) = paste0(colnames(perturations), "_perturb")
combined_perturbations = cbind(normal, perturations)
dim(combined_perturbations)
rownames(metadata_perturb) = colnames(combined_perturbations)
combined_perturbations_seurat = CreateSeuratObject(counts = combined_perturbations, assay = "RNA")
combined_perturbations_seurat[["RNA"]] = CreateAssay5Object(data = combined_perturbations)
combined_perturbations_seurat = AddMetaData(combined_perturbations_seurat, metadata = metadata_perturb)
table(combined_perturbations_seurat$Perturbed)
Idents(combined_perturbations_seurat) = "Perturbed"
perturbed_degs = FindMarkers(combined_perturbations_seurat, only.pos = F, ident.1 = "Perturb",
                             min.pct = 0, min.diff.pct = 0, logfc.threshold = 0)
saveRDS(perturbed_degs, "DEGs/AML_FOSB_Perturbation.rds")


########################################################3
#####CEBPA
########################################################
perturbation_data <- read.table("SCENICPlus_Pipeline/AML/Perturbation/AML_CEBPA_perturbation_iter_11_processed.tsv", header = TRUE, row.names = 1, sep = "\t")
perturbation_data = t(perturbation_data)
aligned_perturbation <- perturbation_data[, colnames(AML_metacells)]
perturbation_matrix <- as.matrix(aligned_perturbation)
perturbation_matrix <- pmax(perturbation_matrix, 0)
AML_metacells[["CEBPA_Perturbation"]] <- CreateAssay5Object(data = perturbation_matrix)
print(AML_metacells)

perturations = as.matrix(GetAssayData(AML_metacells, assay = "CEBPA_Perturbation", layer = "data"))
metadata_perturb = data.frame(Perturbed = c(rep("NoPerturb", dim(AML_metacells)[2]), rep("Perturb", dim(AML_metacells)[2])))
colnames(perturations) = paste0(colnames(perturations), "_perturb")
combined_perturbations = cbind(normal, perturations)
dim(combined_perturbations)
rownames(metadata_perturb) = colnames(combined_perturbations)
combined_perturbations_seurat = CreateSeuratObject(counts = combined_perturbations, assay = "RNA")
combined_perturbations_seurat[["RNA"]] = CreateAssay5Object(data = combined_perturbations)
combined_perturbations_seurat = AddMetaData(combined_perturbations_seurat, metadata = metadata_perturb)
table(combined_perturbations_seurat$Perturbed)
Idents(combined_perturbations_seurat) = "Perturbed"
perturbed_degs = FindMarkers(combined_perturbations_seurat, only.pos = F, ident.1 = "Perturb",
                             min.pct = 0, min.diff.pct = 0, logfc.threshold = 0)
saveRDS(perturbed_degs, "DEGs/AML_CEBPA_Perturbation.rds")


###############################################
DefaultAssay(AML_metacells) = "RNA"
AML_metacells = AddModuleScore(AML_metacells, features = list(hspc_genes), name = "HSPC_TRN_Default")

DefaultAssay(AML_metacells) = "JUNB_Perturbation"
AML_metacells = AddModuleScore(AML_metacells, features = list(hspc_genes), name = "HSPC_TRN_JUNB")

DefaultAssay(AML_metacells) = "HOXA9_Perturbation"
AML_metacells = AddModuleScore(AML_metacells, features = list(hspc_genes), name = "HSPC_TRN_HOXA9")

DefaultAssay(AML_metacells) = "HOXA3_Perturbation"
AML_metacells = AddModuleScore(AML_metacells, features = list(hspc_genes), name = "HSPC_TRN_HOXA3")

DefaultAssay(AML_metacells) = "HOXA5_Perturbation"
AML_metacells = AddModuleScore(AML_metacells, features = list(hspc_genes), name = "HSPC_TRN_HOXA5")

DefaultAssay(AML_metacells) = "FOSB_Perturbation"
AML_metacells = AddModuleScore(AML_metacells, features = list(hspc_genes), name = "HSPC_TRN_FOSB")

DefaultAssay(AML_metacells) = "CEBPA_Perturbation"
AML_metacells = AddModuleScore(AML_metacells, features = list(hspc_genes), name = "HSPC_TRN_CEBPA")

AML_metacells$All = "All"

wilcox.test(AML_metacells$HSPC_TRN_Default1, AML_metacells$HSPC_TRN_JUNB1, alternative = "greater")
wilcox.test(AML_metacells$HSPC_TRN_Default1, AML_metacells$HSPC_TRN_FOSB1, alternative = "greater")
wilcox.test(AML_metacells$HSPC_TRN_Default1, AML_metacells$HSPC_TRN_HOXA31, alternative = "greater")
wilcox.test(AML_metacells$HSPC_TRN_Default1, AML_metacells$HSPC_TRN_HOXA51, alternative = "greater")
wilcox.test(AML_metacells$HSPC_TRN_Default1, AML_metacells$HSPC_TRN_HOXA91, alternative = "greater")
wilcox.test(AML_metacells$HSPC_TRN_Default1, AML_metacells$HSPC_TRN_CEBPA1, alternative = "greater")

p1 = VlnPlot(AML_metacells, features = c("HSPC_TRN_Default1"), group.by = "All", pt.size = 0, same.y.lims = T) + coord_cartesian(ylim = c(0.25, 1.25)) +
  geom_boxplot(width=0.2, fill="white", outlier.shape = NA) + NoLegend() + NoAxes()
p2 = VlnPlot(AML_metacells, features = c("HSPC_TRN_JUNB1"), group.by = "All", pt.size = 0, same.y.lims = T) + coord_cartesian(ylim = c(0.25, 1.25)) +
  geom_boxplot(width=0.2, fill="white", outlier.shape = NA) + NoLegend() + NoAxes()
grid.arrange(p1, p2, nrow = 1)

p1 = VlnPlot(AML_metacells, features = c("HSPC_TRN_Default1"), group.by = "predicted.celltypegeneral_refmap", pt.size = 0, same.y.lims = T) + 
  geom_boxplot(width=0.2, fill="white", outlier.shape = NA)
p2 = VlnPlot(AML_metacells, features = c("HSPC_TRN_JUNB1"), group.by = "predicted.celltypegeneral_refmap", pt.size = 0, same.y.lims = T) + 
  geom_boxplot(width=0.2, fill="white", outlier.shape = NA)
grid.arrange(p1, p2, nrow = 1)

