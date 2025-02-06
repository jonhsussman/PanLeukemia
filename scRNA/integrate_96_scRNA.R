source('scDataAnalysis_Utilities.R')
source('dsAnalysis_utilities.R')
library(Seurat)
library(Signac)
library(ggplot2)
library(parallel)
library(tidyr)
library(dplyr)
library(plyr)
library(readxl)
library(data.table)
library(chromVAR)
library(RColorBrewer)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggprism)
library(patchwork)
library(magrittr)
library(viridis)
library(SCENT)


# Load and Combine all leukemia samples ---------------

files <- dir('PanLeukemia/scRNA/SeuratObj/BySamples/')

files <- files[!(grepl(files, pattern = 'Combined'))]

tmp <- readRDS(paste0('PanLeukemia/scRNA/SeuratObj/BySamples/', files[1]))
mtx <- tmp@assays$RNA@counts
samples = tmp$sample.name
perc.mito = tmp$percent.mt

for(file0 in files[-1]){
  tmp <- readRDS(paste0('PanLeukemia/scRNA/SeuratObj/BySamples/', file0))
  mtx0 <- tmp@assays$RNA@counts
  samples0 = tmp$sample.name
  perc.mito0 = tmp$percent.mt
  
  samples = c(samples, samples0)
  perc.mito = c(perc.mito, perc.mito0)
  
  mtx = t(Matrix.utils::rBind.fill(t(mtx), t(mtx0), fill = 0))
}

rm(tmp, mtx0)


seurat.rna <- CreateSeuratObject(mtx)
seurat.rna <- NormalizeData(seurat.rna)
seurat.rna <- AddMetaData(seurat.rna, metadata = meta)
seurat.rna <- FindVariableFeatures(seurat.rna, nfeatures = 10000)
cell.express.perc = rowSums(seurat.rna@assays$RNA@counts >= 1) /ncol(seurat.rna@assays$RNA@counts)
table(cell.express.perc > 0.01)

sel.features <- setdiff(VariableFeatures(seurat.rna), names(cell.express.perc)[cell.express.perc < 0.01])

VariableFeatures(seurat.rna) = sel.features

cycle3 = fread('PanLeukemia/Signature_GeneList/regev_lab_cell_cycle_genes.txt', header = F)$V1
s.genes = cycle3[1:43]
g2m.genes = cycle3[44:97]
seurat.rna = CellCycleScoring(object = seurat.rna, s.features = s.genes,
                              g2m.features = g2m.genes, set.ident = TRUE)

#heat shock regression
heat_shock_gene = fread('PanLeukemia/Signature_GeneList/heat_shock_geneList.txt')
heat_shock_gene = heat_shock_gene$`Approved symbol`
seurat.rna = AddModuleScore(seurat.rna, features = list(heat_shock_gene),
                            name = 'HeatShock.Score')


seurat.rna <- ScaleData(seurat.rna, features = VariableFeatures(seurat.rna), 
                        vars.to.regress = c("S.Score", "G2M.Score", "perc.mito","HeatShock.Score1", "nCount_RNA"))


seurat.rna <- RunPCA(seurat.rna, npc = 100, verbose = F, 
                     features = VariableFeatures(seurat.rna))

seurat.rna <- RunUMAP(seurat.rna, dims = 1:50, reduction.name = 'UMAP_pca50', reduction.key = 'UMAP_pca50_')

seurat.rna <- FindNeighbors(seurat.rna, dims = 1:50, reduction = 'pca')
seurat.rna <- FindClusters(seurat.rna, res = 0.4)

## filtering variable genes #####

freq_gene <- rowMeans(mtx > 0)
filter.genes = names(which(freq_gene < 0.01))

## reselect variable genes ####
niter = 2
k = 0
topn = 3000
npc = 50
repeat{
  k = k + 1
  if(k > niter) break
  clusters = as.character(seurat.rna$seurat_clusters)
  mtx_by_cls <- sapply(unique(clusters), function(x) {
    
    cl_data <- mtx[, clusters == x]
    
    Matrix::rowSums(cl_data)
    
  })
  mtx_by_cls.norm <- edgeR::cpm(mtx_by_cls, log = T, prior.count = 1)
  sds = sapply(1:nrow(mtx_by_cls.norm), function(x) sd(mtx_by_cls.norm[x, ]))
  names(sds) = rownames(mtx_by_cls.norm)
  sds = sort(sds, decreasing = T)
  sele.genes = names(sds[1:topn])
  sele.genes = setdiff(sele.genes, filter.genes)
  
  VariableFeatures(seurat.rna) <- sele.genes
  
  seurat.rna <- ScaleData(seurat.rna, features = VariableFeatures(seurat.rna), 
                          vars.to.regress =c('S.Score', 'G2M.Score', 'perc.mito',
                                             'HeatShock.Score1', 'nCount_RNA' ))
  seurat.rna <- RunPCA(seurat.rna, npcs = npc, verbose = F, 
                       features = VariableFeatures(seurat.rna))
  
  seurat.rna <- RunUMAP(seurat.rna, dims = 1:npc)
  #seurat.rna <- RunTSNE(seurat.rna, dims = 1:npc)
  seurat.rna <- FindNeighbors(seurat.rna, dims = 1:npc, reduction = 'pca')
  seurat.rna <- FindClusters(seurat.rna, res = 0.4)
  
}
DimPlot(seurat.rna, group.by = 'sample', label = T)

saveRDS(seurat.rna, file = "PanLeukemia/SeuratObj/Combined/seurat_allleukemia_regrCycleHeatShockGenes_pool_VEG3000.rds")


