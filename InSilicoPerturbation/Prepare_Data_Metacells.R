library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(anndata)
library(Signac)
library(Matrix)
library(ggplot2)
library(Signac)
library(data.table)

setwd("/mnt/isilon/tan_lab/sussmanj/Temp/PanLeukemia/SCENICPlus")

#################################
#AML 
#################################
file_ids <- c("143", "150", "209", "2101", "2142", "2144", "2167", "31", "1832")
sample_id <- "7767"
seurat_objects <- list()
for (id in file_ids) {
  print(id)
  rna_path <- paste0("../Coembed_Results_Metacell_Counts/rna_metacell_mtx_", sample_id, "-", id, "_k.anchor30.rds")
  atac_path <- paste0("../Coembed_Results_Metacell_Counts/atac_metacell_mtx_", sample_id, "-", id, "_k.anchor30.rds")
  rna_mtx <- readRDS(rna_path)
  atac_mtx <- readRDS(atac_path)
  seurat_obj <- CreateSeuratObject(counts = rna_mtx, assay = "RNA")
  seurat_obj[["ATAC"]] <- CreateChromatinAssay(counts = atac_mtx)
  seurat_obj$orig.ident <- paste0(sample_id, "-", id)
  seurat_objects[[paste0("seurat_", sample_id, "_", id)]] <- seurat_obj
}

AML_Seurat_metacells = merge(seurat_objects[[1]], seurat_objects[-1])
DefaultAssay(AML_Seurat_metacells) = "RNA"
AML_Seurat_metacells = JoinLayers(AML_Seurat_metacells)

#Unintegrated analysis
AML_Seurat_metacells <- NormalizeData(AML_Seurat_metacells)
AML_Seurat_metacells <- FindVariableFeatures(AML_Seurat_metacells, nfeatures = 2000)
AML_Seurat_metacells <- ScaleData(AML_Seurat_metacells) 
AML_Seurat_metacells <- RunPCA(AML_Seurat_metacells, verbose = T) 
AML_Seurat_metacells <- RunUMAP(AML_Seurat_metacells, reduction = 'pca', dim = 1:30)
DimPlot(AML_Seurat_metacells, group.by = "orig.ident") + coord_fixed()

#Project to reference 
healthy_ref = readRDS('/mnt/isilon/tan_lab/sussmanj/Temp/T_ALL/Healthy_Reference_Fixed.rds')
ob.list <- seurat_objects
for (i in 1:length(ob.list)) {
  print(i)
  DefaultAssay(ob.list[[i]]) = "RNA"
  ob.list[[i]] <- NormalizeData(ob.list[[i]])
  ob.list[[i]] <- FindVariableFeatures(ob.list[[i]])
  ob.list[[i]] <- ScaleData(ob.list[[i]])
  ob.list[[i]] <- RunPCA(ob.list[[i]], npcs = 30)
  anchors <- FindTransferAnchors(
    reference = healthy_ref,
    query = ob.list[[i]],
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    dims = 1:30
  )
  ob.list[[i]] <- MapQuery(
    anchorset = anchors, 
    query = ob.list[[i]],
    reference = healthy_ref,
    refdata = list(
      t.traj.ptime = "curve2", 
      m.traj.ptime = "curve5",
      celltype_refmap = "cell.type", 
      celltypegeneral_refmap = "cell.type.general", 
      trajectory = "trajectory"
    ),
    reference.reduction = "pca", 
    reduction.model = "umap.model", 
    transferdata.args = list(k.weight = 40)
  )
}
projected_seurat = merge(ob.list[[1]], y=ob.list[-1], merge.dr = c("ref.umap"))
DimPlot(projected_seurat, reduction = "ref.umap", group.by = "predicted.celltype_refmap", label = T) + coord_fixed()
saveRDS(projected_seurat, "AML_Metacells_Projected.rds")

#################################
#T-ALL
#################################
file_ids <- c(
  "PASKMG", "PASZKM", "PASZMC", "PATIPB", "PATMYZ", "PATTDP", "PAUFAM", "PAUHWY",
  "PAUNDK", "PAURIX", "PAURXZ", "PAUUZY", "PAVFKN", "PAVINC", "PAVJCH", "PAVLII",
  "PAVLKA", "PAVLUN", "PAVSEI", "PAVSRU", "PAVTCV", "PAVTXP", "PAVVVF", "PAVYVY",
  "PAWGWD", "PASWWT", "PATENL", "PATEVG", "PATPKZ", "PAUKIZ", "PATDFE", "PAUMAV",
  "PAUMXB", "PAUNZE", "PAUPTX", "PAUPVR", "PAUYJE", "PAVUFK", "PAVVVK", "PAWIIR"
)
seurat_objects <- list()
for (id in file_ids) {
  print(id)
  rna_path <- paste0("../Coembed_Results_Metacell_Counts/rna_metacell_mtx_", id, "_k.anchor30.rds")
  atac_path <- paste0("../Coembed_Results_Metacell_Counts/atac_metacell_mtx_", id, "_k.anchor30.rds")
  rna_mtx <- readRDS(rna_path)
  atac_mtx <- readRDS(atac_path)
  seurat_obj <- CreateSeuratObject(counts = rna_mtx, assay = "RNA")
  seurat_obj[["ATAC"]] <- CreateChromatinAssay(counts = atac_mtx)
  seurat_obj$orig.ident <- paste0(sample_id, "-", id)
  seurat_objects[[paste0("seurat_", sample_id, "_", id)]] <- seurat_obj
}

TALL_Seurat_metacells = merge(seurat_objects[[1]], seurat_objects[-1])
DefaultAssay(TALL_Seurat_metacells) = "RNA"
TALL_Seurat_metacells = JoinLayers(TALL_Seurat_metacells)

#Unintegrated analysis
TALL_Seurat_metacells <- NormalizeData(TALL_Seurat_metacells)
TALL_Seurat_metacells <- FindVariableFeatures(TALL_Seurat_metacells, nfeatures = 2000)
TALL_Seurat_metacells <- ScaleData(TALL_Seurat_metacells) 
TALL_Seurat_metacells <- RunPCA(TALL_Seurat_metacells, verbose = T) 
TALL_Seurat_metacells <- RunUMAP(TALL_Seurat_metacells, reduction = 'pca', dim = 1:30)
DimPlot(TALL_Seurat_metacells, group.by = "orig.ident") + coord_fixed()

#Project to reference 
healthy_ref = readRDS('/mnt/isilon/tan_lab/sussmanj/Temp/T_ALL/Healthy_Reference_Fixed.rds')
ob.list <- seurat_objects
ob.list[["seurat_7767_PAUYJE"]] = NULL
for (i in 1:length(ob.list)) {
  print(i)
  DefaultAssay(ob.list[[i]]) = "RNA"
  ob.list[[i]] <- NormalizeData(ob.list[[i]])
  ob.list[[i]] <- FindVariableFeatures(ob.list[[i]])
  ob.list[[i]] <- ScaleData(ob.list[[i]])
  ob.list[[i]] <- RunPCA(ob.list[[i]], npcs = 30)
  anchors <- FindTransferAnchors(
    reference = healthy_ref,
    query = ob.list[[i]],
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    dims = 1:30
  )
  ob.list[[i]] <- MapQuery(
    anchorset = anchors, 
    query = ob.list[[i]],
    reference = healthy_ref,
    refdata = list(
      t.traj.ptime = "curve2", 
      m.traj.ptime = "curve5",
      celltype_refmap = "cell.type", 
      celltypegeneral_refmap = "cell.type.general", 
      trajectory = "trajectory"
    ),
    reference.reduction = "pca", 
    reduction.model = "umap.model", 
    transferdata.args = list(k.weight = 40)
  )
}
projected_seurat = merge(ob.list[[1]], y=ob.list[-1], merge.dr = c("ref.umap"))
DimPlot(projected_seurat, reduction = "ref.umap", group.by = "predicted.celltype_refmap", label = T) + coord_fixed()
saveRDS(projected_seurat, "TALL_Metacells_Projected.rds")

#################################
#B-ALL 
#################################
for(i in 1:length(seurat_objects)){
  print(head(rownames(seurat_objects[[i]]@assays$ATAC)))
}
    
file_ids <- c(
  "1154", "1979", "2184", "2456", "2744", "3315", "PAYKGI", "PAYLNH", "PAYSBA", 
  "PAYUZJ", "PAYUZM", "PAYWJZ", "PAYWKL", "PAYYBG", "PAYYNY", "PAYZLC", "PAYZVY", 
  "PAYZWN", "PAZBGV", "PAZBLA", "PAZBSZ", "PAZFPH", "PAZGKI", "PAPIWL", "PAREZH", 
  "PATGLR", "PATSBP", "PAVALZ", "PAVXLK", "PAXEIP", "PAXFFK", "PAXYRK", "PAYFXM", 
  "PAYLZC", "PAYPDR"
)
seurat_objects <- list()
for (id in file_ids) {
  print(id)
  rna_path <- paste0("../Coembed_Results_Metacell_Counts/rna_metacell_mtx_", id, "_k.anchor30.rds")
  atac_path <- paste0("../Coembed_Results_Metacell_Counts/atac_metacell_mtx_", id, "_k.anchor30.rds")
  rna_mtx <- readRDS(rna_path)
  atac_mtx <- readRDS(atac_path)
  seurat_obj <- CreateSeuratObject(counts = rna_mtx, assay = "RNA")
  seurat_obj[["ATAC"]] <- CreateChromatinAssay(counts = atac_mtx)
  seurat_obj$orig.ident <- paste0(id)
  seurat_objects[[paste0("seurat_", id)]] <- seurat_obj
}

BALL_Seurat_metacells = merge(seurat_objects[[1]], seurat_objects[-1])
DefaultAssay(BALL_Seurat_metacells) = "RNA"
BALL_Seurat_metacells = JoinLayers(BALL_Seurat_metacells)

#Unintegrated analysis
BALL_Seurat_metacells <- NormalizeData(BALL_Seurat_metacells)
BALL_Seurat_metacells <- FindVariableFeatures(BALL_Seurat_metacells, nfeatures = 2000)
BALL_Seurat_metacells <- ScaleData(BALL_Seurat_metacells) 
BALL_Seurat_metacells <- RunPCA(BALL_Seurat_metacells, verbose = T) 
BALL_Seurat_metacells <- RunUMAP(BALL_Seurat_metacells, reduction = 'pca', dim = 1:30)
DimPlot(BALL_Seurat_metacells, group.by = "orig.ident") + coord_fixed()

#Project to reference 
healthy_ref = readRDS('/mnt/isilon/tan_lab/sussmanj/Temp/T_ALL/Healthy_Reference_Fixed.rds')
ob.list <- seurat_objects
ob.list[["seurat_PAYZVY"]] = NULL

for (i in 1:length(ob.list)) {
  print(i)
  DefaultAssay(ob.list[[i]]) = "RNA"
  ob.list[[i]] <- NormalizeData(ob.list[[i]])
  ob.list[[i]] <- FindVariableFeatures(ob.list[[i]])
  ob.list[[i]] <- ScaleData(ob.list[[i]])
  ob.list[[i]] <- RunPCA(ob.list[[i]], npcs = 30)
  anchors <- FindTransferAnchors(
    reference = healthy_ref,
    query = ob.list[[i]],
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    dims = 1:30
  )
  ob.list[[i]] <- MapQuery(
    anchorset = anchors, 
    query = ob.list[[i]],
    reference = healthy_ref,
    refdata = list(
      t.traj.ptime = "curve2", 
      m.traj.ptime = "curve5",
      celltype_refmap = "cell.type", 
      celltypegeneral_refmap = "cell.type.general", 
      trajectory = "trajectory"
    ),
    reference.reduction = "pca", 
    reduction.model = "umap.model", 
    transferdata.args = list(k.weight = 40)
  )
}
projected_seurat = merge(ob.list[[1]], y=ob.list[-1], merge.dr = c("ref.umap"))
DimPlot(projected_seurat, reduction = "ref.umap", group.by = "predicted.celltype_refmap", label = T) + coord_fixed()
saveRDS(projected_seurat, "BALL_Metacells_Projected.rds")


####################
#Save data for SCENIC+
####################
#Save regions to bed format 
convert_to_bed <- function(region) {
  parts <- unlist(strsplit(region, "-"))
  chr <- parts[1]
  start <- as.integer(parts[2])
  end <- as.integer(parts[3])
  bed <- c(chr, start, end)
  return(bed)
}

#####AML
aml_metacells_seurat = readRDS("AML_Metacells_Projected.rds")

#Save RNA-seq
DefaultAssay(aml_metacells_seurat) <- "RNA"
aml_metacells_seurat = JoinLayers(aml_metacells_seurat)
aml.rna.loom <- as.loom(aml_metacells_seurat, filename = "Data_SCENICplus/AML_Metacells_RNA.loom", verbose = TRUE) 
aml.rna.loom$close_all() 

#Save ATAC-seq data 
counts_matrix <- as.matrix(aml_metacells_seurat@assays$ATAC@counts)
counts_sparse <- Matrix::Matrix(counts_matrix , sparse = T)
writeMM(counts_sparse, file = "Data_SCENICplus/AML_Metacells_ATAC_Peaks_Sparse.mtx")
cell_names <- colnames(counts_sparse)
region_names <- rownames(counts_sparse)
cell_names_file <- "Data_SCENICplus/AML_Metacells_ATAC_Cell_Names.txt"
region_names_file <- "Data_SCENICplus/AML_Metacells_ATAC_Region_Names.txt"
metadata_frame <- aml_metacells_seurat@meta.data
write.table(cell_names, file = cell_names_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(region_names, file = region_names_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(metadata_frame, file = "Data_SCENICplus/AML_Metacells_ATAC_Metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)
bed_data <- apply(as.matrix(region_names), 1, convert_to_bed)
write.table(t(bed_data), "Data_SCENICplus/AML_Metacells_ATAC_Region_Names.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#####T-ALL
tall_metacells_seurat = readRDS("TALL_Metacells_Projected.rds")
tall_metacells_seurat = JoinLayers(tall_metacells_seurat)

#Save RNA-seq
DefaultAssay(tall_metacells_seurat) <- "RNA"
tall.rna.loom <- as.loom(tall_metacells_seurat, filename = "Data_SCENICplus/TALL_Metacells_RNA.loom", verbose = TRUE) 
tall.rna.loom$close_all() 

#Save ATAC-seq data 
counts_matrix <- as.matrix(tall_metacells_seurat@assays$ATAC@counts)
counts_sparse <- Matrix::Matrix(counts_matrix , sparse = T)
writeMM(counts_sparse, file = "Data_SCENICplus/TALL_Metacells_ATAC_Peaks_Sparse.mtx")
cell_names <- colnames(counts_sparse)
region_names <- rownames(counts_sparse)
cell_names_file <- "Data_SCENICplus/TALL_Metacells_ATAC_Cell_Names.txt"
region_names_file <- "Data_SCENICplus/TALL_Metacells_ATAC_Region_Names.txt"
metadata_frame <- tall_metacells_seurat@meta.data
write.table(cell_names, file = cell_names_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(region_names, file = region_names_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(metadata_frame, file = "Data_SCENICplus/TALL_Metacells_ATAC_Metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)
bed_data <- apply(as.matrix(region_names), 1, convert_to_bed)
write.table(t(bed_data), "Data_SCENICplus/TALL_Metacells_ATAC_Region_Names.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#####B-ALL
ball_metacells_seurat = readRDS("BALL_Metacells_Projected.rds")
ball_metacells_seurat = JoinLayers(ball_metacells_seurat)

#Save RNA-seq
DefaultAssay(ball_metacells_seurat) <- "RNA"
ball.rna.loom <- as.loom(ball_metacells_seurat, filename = "Data_SCENICplus/BALL_Metacells_RNA.loom", verbose = TRUE) 
ball.rna.loom$close_all() 

#Save ATAC-seq data 
counts_matrix <- as.matrix(ball_metacells_seurat@assays$ATAC@counts)
counts_sparse <- Matrix::Matrix(counts_matrix , sparse = T)
writeMM(counts_sparse, file = "Data_SCENICplus/BALL_Metacells_ATAC_Peaks_Sparse.mtx")
cell_names <- colnames(counts_sparse)
region_names <- rownames(counts_sparse)
cell_names_file <- "Data_SCENICplus/BALL_Metacells_ATAC_Cell_Names.txt"
region_names_file <- "Data_SCENICplus/BALL_Metacells_ATAC_Region_Names.txt"
metadata_frame <- ball_metacells_seurat@meta.data
write.table(cell_names, file = cell_names_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(region_names, file = region_names_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(metadata_frame, file = "Data_SCENICplus/BALL_Metacells_ATAC_Metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)
bed_data <- apply(as.matrix(region_names), 1, convert_to_bed)
write.table(t(bed_data), "Data_SCENICplus/BALL_Metacells_ATAC_Region_Names.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

