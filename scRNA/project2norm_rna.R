source('scDataAnalysis_Utilities.R')

## build the normal tranjectory (in process_healthydonor_final.R) ####
seurat.rna <- readRDS('Seurat_Objects/scRNA/seurat_pool_logNorm_gini_FiveHD_10Xv3_downsample10000HSPC.rds')


## revovery all the previous transformation done for thess healthy donor cells ####
vegs = VariableFeatures(seurat.rna)
feature.load = seurat.rna@reductions$pca@feature.loadings
# get the coefficient from ScaleData for each gene
coef.lm = centers.vegs = sd.vegs = list()
data0 = seurat.rna@assays$RNA@data
mdata0 = seurat.rna@meta.data
for(gene0 in vegs){
  expr0 = data0[gene0, ]
  perc.mito0 = mdata0[, 'perc.mito']
  numi0 = mdata0[, 'nCount_RNA']
  model0 = lm(expr0 ~ perc.mito0 + numi0)
  residual0 = residuals(model0)
  coef.lm[[gene0]] = coef(model0)
  
  centers.vegs[[gene0]] = mean(residual0)
  sd.vegs[[gene0]] = sd(residual0)
}
coef.lm = do.call('rbind', coef.lm)
centers.vegs = do.call('c', centers.vegs)
sd.vegs = do.call('c', sd.vegs)

## save coefficients for future use
save(coef.lm, centers.vegs, sd.vegs, file = 
       'Seurat_Objects/scRNA/coef_hd.RData')

## recovery the umap model
set.seed(42) 
npc = 20
umap_model = uwot::umap(as.matrix(seurat.rna@reductions$pca@cell.embeddings[, 1:npc]),
                        n_neighbors = 30, min_dist = 0.3,
                        ret_model = T, metric = 'cosine')

health_ann = seurat.rna@meta.data
health_ann = subset(health_ann, select = c('sample', 'nCount_RNA', 'nFeature_RNA', 'perc.mito',  'Ctype'))
health_ann = data.table(health_ann, keep.rownames = T)
setkey(health_ann, rn)


## direct projection ####
## normalize the patient data the same way as did for healthy donors
seurat.all.patients = readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')
mllr_ann = data.table(seurat.all.patients@meta.data, keep.rownames = T)
mllr_ann$Ctype = mllr_ann$Ctype0

cell.map = NULL
for(sampleName in unique(mllr_ann$sample)){
  seurat.mllr = subset(seurat.all.patients, sample == sampleName)
  mllr_ann1 = subset(mllr_ann, sample == sampleName)
  setkey(mllr_ann1, rn)
  
  seurat.mllr = NormalizeData(seurat.mllr)
  
  ## scale the data using coefficent learned in healthy donor cells
  shared.genes = intersect(rownames(seurat.mllr), vegs)
  coef.lm1 = coef.lm[shared.genes, ]
  centers.vegs1 = centers.vegs[shared.genes]
  sd.vegs1 = sd.vegs[shared.genes]
  
  
  mtx1 = seurat.mllr@assays$RNA@data[shared.genes, ]
  mdata1 = seurat.mllr@meta.data
  perc.mito1 = mdata1[, 'perc.mito']
  nCount_RNA1 = mdata1[, 'nCount_RNA']
  covMat = cbind(rep(1, length(perc.mito1)), perc.mito1, nCount_RNA1)
  mtx1.scale = (mtx1 - coef.lm1 %*% t(covMat) - centers.vegs1)/sd.vegs1
  
  
  U = feature.load[shared.genes, 1:npc]
  
  pca1 = t(mtx1.scale) %*% U ##projected pcs
  
  
  ## to be consistent with seurat seed
  umap1 = uwot::umap_transform(as.matrix(pca1), umap_model)
  
  rownames(umap1) = rownames(pca1)
  umap1 = umap1[rownames(umap1) %in% mllr_ann1$rn, ]
  Ctype1 = mllr_ann1[rownames(umap1), ]$Ctype
  Ctype1 = sapply(Ctype1, function(x) gsub('_Cells', '', x))
  
  pdata = rbind(umap_model$embedding, umap1)
  pdata = data.table(pdata)
  pdata = cbind(pdata, c(rep('Healthy', nrow(umap_model$embedding)), 
                         Ctype1))
  names(pdata) = c('UMAP_1', 'UMAP_2', 'ctype')
  
  myColors = brewer.pal(8, name = 'Set1')
  names(myColors) <- c('Blasts', 'Monocytes', 'Mature_B' , 'Progenitors', 
                       'NKT', 'pDC',  'T',  'Healthy')
  myColors['Healthy'] = '#cccccc'
  
  pdata$ctype = factor(pdata$ctype, levels = c('Healthy', sort(unique(Ctype1))) )
  pp1 <- ggplot(pdata, aes(x = UMAP_1, y = UMAP_2, col = ctype)) + 
    geom_point(size = 0.2) + 
    scale_colour_manual(name = "annotation", values = myColors)
  
  ggsave(pp1, filename = paste0('Figures/scRNA/proj2Healthy/final/proj_direct_', sampleName,
                                '.eps'), device = 'eps', width = 9, height = 7)
  
  
  ## assign the cell type from nearest healthy donor cells
  rownames(umap_model$embedding) = colnames(seurat.rna)
  
  ss.sub = pracma::distmat(umap1, umap_model$embedding)
  
  
  mllr2health <- sapply(1:nrow(umap1), function(x) names(which.min(ss.sub[x, ])))
  
  cell.map.tmp <- data.table('patient_bc' = rownames(umap1),
                             'healthy_bc' = mllr2health)
  cell.map = rbind(cell.map, cell.map.tmp)
  
}


map2fname = 'MetaData/scRNA/MLLrCells_proj2healthyCells_final.csv'
setkey(mllr_ann, rn)

cell.map[, 'ctype_healthy' := health_ann[J(cell.map$healthy_bc)]$Ctype]
cell.map[, 'sample_healthy' := health_ann[J(cell.map$healthy_bc)]$sample]
cell.map[, 'ctype_patient' := mllr_ann[J(cell.map$patient_bc)]$Ctype]
cell.map[, 'sample_patient' := mllr_ann[J(cell.map$patient_bc)]$sample]

write.table(cell.map, file = map2fname, sep = '\t',
            row.names = F, quote = F)

## add projected cell type information to patient seurat
proj_ctype = subset(cell.map, select = c('ctype_healthy', 'healthy_bc', 'ctype_patient'))
proj_ctype$projCtype1 <- proj_ctype$ctype_healthy
proj_ctype$projCtype_merge1 <- proj_ctype$ctype_healthy
proj_ctype$projCtype_merge1 <- ifelse(proj_ctype$projCtype_merge1 %in% c('HSPC', 'LMPP'),
                                     'Early-prog', proj_ctype$projCtype_merge1)
proj_ctype$projCtype_merge1 <- ifelse(proj_ctype$projCtype_merge1 %in% c('Pre-pro-B', 'Pro-B'),
                                     'Pro-B', proj_ctype$projCtype_merge1)

proj_ctype$projCtype1 = paste0(proj_ctype$projCtype1, '-like')
proj_ctype$projCtype_merge1 = paste0(proj_ctype$projCtype_merge1, '-like')

proj_ctype = data.frame(proj_ctype)
rownames(proj_ctype) = cell.map$patient_bc
names(proj_ctype)[2] = 'healthy_bc1'
seurat.all.patients <- AddMetaData(seurat.all.patients, 
                                   metadata = proj_ctype)
saveRDS(seurat.all.patients, file = 'Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')
