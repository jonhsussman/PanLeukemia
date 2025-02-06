## process data from healty donors 

source('scDataAnalysis_Utilities.R')


## gini selection of VEG ####
files = dir('PanLeukemia/scRNA/Seurat_Objects/HealthyDonor/')

tmp <- readRDS(paste0('PanLeukemia/scRNA/Seurat_Objects/HealthyDonor/', files[1]))
tmp <- subset(tmp, nCount_RNA < 40000)
mtx <- tmp@assays$RNA@counts
samples = tmp$sample
perc.mito = tmp$perc.mito

for(file0 in files[-1]){
  tmp <- readRDS(paste0('Seurat_Objects/HealthyDonor/', file0))
  tmp <- subset(tmp, nCount_RNA < 40000 )
  mtx0 <- tmp@assays$RNA@counts
  samples0 = tmp$sample
  perc.mito0 = tmp$perc.mito
  
  samples = c(samples, samples0)
  perc.mito = c(perc.mito, perc.mito0)
  
  mtx = t(Matrix.utils::rBind.fill(t(mtx), t(mtx0), fill = 0))
}

rm(tmp, mtx0)

ids_cd34p = grep(samples, pattern = '_CD34')
ids_cd38n = grep(samples, pattern = 'CD38')

set.seed(2019)
ids_sub = sample(c(ids_cd34p, ids_cd38n), 10000)
ids_other = setdiff(1:length(samples), c(ids_cd34p, ids_cd38n))
ids_sele = sort(c(ids_other, ids_sub))
mtx_sub = mtx[, ids_sele]

seurat.rna <- CreateSeuratObject(mtx_sub)
seurat.rna <- NormalizeData(seurat.rna)
seurat.rna$sample = samples[ids_sele]
seurat.rna$perc.mito = perc.mito[ids_sele]
seurat.rna <- FindVariableFeatures(seurat.rna, nfeatures = 1000)

seurat.rna <- ScaleData(seurat.rna, features = VariableFeatures(seurat.rna), 
                        vars.to.regress = c('perc.mito', 'nCount_RNA'))
seurat.rna <- RunPCA(seurat.rna, npc = 30, verbose = F, 
                     features = VariableFeatures(seurat.rna))

seurat.rna <- RunUMAP(seurat.rna, dims = 1:20)
#seurat.rna <- RunTSNE(seurat.rna, dims = 1:20)
seurat.rna <- FindNeighbors(seurat.rna, dims = 1:20, reduction = 'pca')
seurat.rna <- FindClusters(seurat.rna, res = 0.1)
DimPlot(seurat.rna, group.by = 'sample')

niter = 2
k = 0
repeat{
  k = k + 1
  if(k > niter) break
  gini_genes <- ifg_select(seurat.rna@assays$RNA@counts, 
                           seurat.rna$seurat_clusters, 
                           gini_cut_qt = 0.9)$include_g
  
  VariableFeatures(seurat.rna) <- gini_genes
  
  seurat.rna <- ScaleData(seurat.rna, features = VariableFeatures(seurat.rna), 
                          vars.to.regress = c('perc.mito', 'nCount_RNA'))
  seurat.rna <- RunPCA(seurat.rna, npc = 30, verbose = F, 
                       features = VariableFeatures(seurat.rna))
  
  seurat.rna <- RunUMAP(seurat.rna, dims = 1:20)
  #seurat.rna <- RunTSNE(seurat.rna, dims = 1:20)
  seurat.rna <- FindNeighbors(seurat.rna, dims = 1:20, reduction = 'pca')
  seurat.rna <- FindClusters(seurat.rna, res = 0.1)
  
}
seurat.rna <- RunPCA(seurat.rna, npc = 100, verbose = F, 
                     features = VariableFeatures(seurat.rna))
seurat.rna <- RunTSNE(seurat.rna, dims = 1:20)
seurat.rna <- RunUMAP(seurat.rna, dims = 1:50)
DimPlot(seurat.rna, group.by = 'sample')


FeaturePlot(seurat.rna, c('CD34', 'CD38', 'CD14', 'FCGR3A', 'CD19', 'MS4A1',
                          'CD3E', 'GNLY', 'GATA1'))

FeaturePlot(seurat.rna, c('IGLL1', 'DNTT'))
FeaturePlot(seurat.rna, c('CD1C', 'IL3RA'))
FeaturePlot(seurat.rna, c('MPO', 'ELANE'))


saveRDS(seurat.rna, 'PanLeukemia/scRNA/Seurat_Objects/HealthyDonor/seurat_pool_logNorm_gini_HDRef_downsample10000HSPC.rds')


