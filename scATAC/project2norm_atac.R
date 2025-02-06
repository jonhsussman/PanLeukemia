source('scDataAnalysis_Utilities.R')

## build the normal tranjectory (in process_healthydonor_scATAC.R) ####
seuratAtacPath = 'Seurat_Objects/scATAC/seurat_9HDsamplesPlusThreeHD_TFIDF_usingDAPs1.rds'

## set base name for automatically naming output files
bname = basename(seuratAtacPath)
bname = gsub('.rds', '', bname)

seurat.hd <- readRDS(seuratAtacPath)

health_ann = seurat.hd@meta.data
health_ann = data.table(health_ann, keep.rownames = T)
setkey(health_ann, rn)


## revovery all the previous transformation done for thess healthy donor cells ####
vaps = VariableFeatures(seurat.hd)
feature.load = seurat.hd@reductions$pca@feature.loadings



# get the tf-idf 
mtx0 = 1*(seurat.hd@assays$ATAC@counts > 0)
mtx0 = mtx0[vaps, ]

npeaks <- colSums(x = mtx0)
tf <- t(t(mtx0) / npeaks)
idf <- ncol(mtx0) / rowSums(x = mtx0)
norm.mtx0 <- Diagonal(n = length(idf), x = idf) %*% tf
all(names(idf) == vaps)

# get the coefficient from ScaleData for each gene
coef.lm = centers.vaps = sd.vaps = list()
data0 = seurat.atac@assays$ATAC@data

mdata0 = seurat.hd@meta.data
for(peak0 in vaps){
  expr0 = data0[peak0, ]
  #numi0 = mdata0$nCount_ATAC
  #model0 = lm(expr0 ~ numi0)
  model0 = lm(expr0 ~ 1)  ## no confounder to regress out
  residual0 = residuals(model0)
  coef.lm[[peak0]] = coef(model0)
  
  centers.vaps[[peak0]] = mean(residual0)
  sd.vaps[[peak0]] = sd(residual0)
}
coef.lm = do.call('rbind', coef.lm)
centers.vaps = do.call('c', centers.vaps)
sd.vaps = do.call('c', sd.vaps)

## save coefficient for future use

save(coef.lm, centers.vaps, sd.vaps, file = 
       paste0('Seurat_Objects/scATAC/coef_', bname, '.RData'))

## recovery the umap model
npc = ncol(seurat.hd@reductions$pca@cell.embeddings)
npc = 20
set.seed(42) 
umap_model = uwot::umap(as.matrix(seurat.hd@reductions$pca@cell.embeddings[, 1:npc]),
                        n_neighbors = 30, min_dist = 0.3,
                        ret_model = T, metric = 'cosine')



## direct projection of Infant MLLr ####
## normalize the patient data the same way as did for healthy donors

nvap = 10000
fileName = paste0('Seurat_Objects/scATAC/seurat_pool_18MLLr_TFIDF_vap', nvap, '.rds')
seurat.mllr = readRDS(fileName)
mllr_ann = seurat.mllr@meta.data

mllr_ann = data.table(mllr_ann, keep.rownames = T)
mllr_ann$rn0 = mllr_ann$rn
mllr_ann$rn = sapply(mllr_ann$rn, function(x){
  len = nchar(x)
  substr(x, len - 15, len)
}) ## remove sample id from cell barcode (to be consistent with colname of filtered mtx)

mllr_ann$bc_sample = paste0(mllr_ann$sample, '_', mllr_ann$rn)

cell.map = NULL
dir0 = '/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/MLLr/'
for(sampleName in unique(mllr_ann$sample)){
  #mtx <- read_mtx_scATACpro(paste0(dir0, sampleName, '/output/matrix4proj/matrix.mtx'))
  ## annotate rownames
  #tss_ann = fread('/mnt/isilon/tan_lab/yuw1/scATAC-pro/annotation/hg38_tss.bed')
  #names(tss_ann)[c(1:4, 6, 7)] = c('chr', 'start', 'end', 'gene_name', 'strand',
  #                                 'gene_type')
  #mtx = assignGene2Peak(mtx, gene_ann = tss_ann)
  ## save to reuse
  #saveRDS(mtx, paste0('matrix/scATAC/matrix4proj_', sampleName, '_batch1WithPubAtac.rds'))
  
  mtx = readRDS(paste0('matrix/scATAC/matrix4proj_', sampleName, '_batch1WithPubAtac.rds'))
  mllr_ann1 = subset(mllr_ann, sample == sampleName)
  setkey(mllr_ann1, rn)
  
  all(colnames(mtx) %in% mllr_ann1$rn)
  
  
  ## scale the data using coefficent learned in healthy donor cells
  shared.peaks = intersect(rownames(mtx), vaps)
  coef.lm1 = coef.lm[shared.peaks, ]
  centers.vaps1 = centers.vaps[shared.peaks]
  sd.vaps1 = sd.vaps[shared.peaks]
  idf1 = idf[shared.peaks]
  
  mtx1 = 1*(mtx[shared.peaks, ] > 0)
  covMat = cbind(rep(1, ncol(mtx1)))
  mtx1 = t(t(mtx1)/colSums(mtx1))  # tf
  mtx1 = Diagonal(n = length(idf1), x = idf1) %*% mtx1
  mtx1.scale = (mtx1 - coef.lm1 %*% t(covMat) - centers.vaps1)/sd.vaps1
  
  
  U = feature.load[shared.peaks, 1:npc]
  
  pca.mllr = t(mtx1.scale) %*% U ##projected pcs
  
  
  ## to be consistent with seurat seed
  umap1 = uwot::umap_transform(as.matrix(pca.mllr), umap_model)
  
  rownames(umap1) = rownames(pca.mllr)
  umap1 = umap1[rownames(umap1) %in% mllr_ann1$rn, ]
  Ctype1 = mllr_ann1[rownames(umap1), ]$Ctype0
  
  pdata = rbind(umap_model$embedding, umap1)
  pdata = data.table(pdata)
  pdata = cbind(pdata, c(rep('Healthy', nrow(umap_model$embedding)), 
                         Ctype1))
  names(pdata) = c('UMAP_1', 'UMAP_2', 'ctype')
  
  myColors = brewer.pal(7, name = 'Set1')
  names(myColors) <- c('Blasts', 'Monocytes', 'Mature_B' , 'Progenitors', 
                       'T/NK', 'Healthy',   'pDC')
  myColors['Healthy'] = '#cccccc'
  
  pdata$ctype = factor(pdata$ctype, levels = c('Healthy', sort(unique(Ctype1))) )
  pp1 <- ggplot(pdata, aes(x = UMAP_1, y = UMAP_2, col = ctype)) + 
    geom_point(size = 0.2) + 
    scale_colour_manual(name = "annotation", values = myColors)
  
  ggsave(pp1, device = 'eps', width = 9, height = 7,
         filename = paste0('Figures/scATAC/proj2Healthy_final/proj_direct_', 
                           sampleName, '.eps'))
  
  
  ## assign the cell type from nearest healthy donor cells
  rownames(umap_model$embedding) = colnames(seurat.hd)
  
  ss.sub = pracma::distmat(umap1, umap_model$embedding)
  #ss.sub = pracma::distmat(as.matrix(pca.mllr), pca.hd)
  
  if(useHighConfLabel) ss.sub = ss.sub[, highconf.bcs]
  
  mllr2health <- sapply(1:nrow(umap1), function(x) names(which.min(ss.sub[x, ])))
  
  cell.map.tmp <- data.table('patient_bc' = rownames(umap1),
                             'healthy_bc' = mllr2health,
                             'sample_patient' = sampleName)
  cell.map = rbind(cell.map, cell.map.tmp)
  
}

cell.map[, 'ctype_healthy' := health_ann[J(cell.map$healthy_bc)]$seurat_ctype]
cell.map[, 'sample_healthy' := health_ann[J(cell.map$healthy_bc)]$sample]
cell.map[, 'bc_sample' := paste0(sample_patient,'_', patient_bc)]
setkey(mllr_ann, bc_sample)
cell.map[, 'ctype_patient' := mllr_ann[J(cell.map$bc_sample)]$Ctype0]

table(cell.map[ctype_patient == 'Progenitors']$ctype_healthy)
table(cell.map[ctype_patient == 'Blasts']$ctype_healthy)

map2fname = paste0('MetaData/scATAC/MLLrCells_proj2healthyCells_', 
                   bname, '.csv')

write.table(cell.map, file = map2fname, sep = '\t',
            row.names = F, quote = F)

## save to patient seurat
setkey(cell.map, bc_sample)
setkey(mllr_ann, rn0)
seurat.mllr$bc_sample = mllr_ann[colnames(seurat.mllr)]$bc_sample

cell.map = cell.map[seurat.mllr$bc_sample, ]
proj_ctype = subset(cell.map, select = c('ctype_healthy', 'healthy_bc', 'ctype_patient'))
proj_ctype$projCtype1 = proj_ctype$ctype_healthy
proj_ctype$projCtype_merge1 <- ifelse(proj_ctype$projCtype1  %in% c('HSPC', 'LMPP'),
                                     'Early-prog', proj_ctype$projCtype1)
proj_ctype$projCtype_merge1 <- ifelse(proj_ctype$projCtype1  %in% c('Pre-pro-B', 'Pro-B'),
                                     'Pro-B', proj_ctype$projCtype_merge1)

proj_ctype$projCtype1 = paste0(proj_ctype$projCtype1, '-like')
proj_ctype$projCtype_merge1 = paste0(proj_ctype$projCtype_merge1, '-like')

proj_ctype = data.frame(proj_ctype)
rownames(proj_ctype) = colnames(seurat.mllr)

seurat.mllr <- AddMetaData(seurat.mllr, 
                           metadata = proj_ctype)

saveRDS(seurat.mllr, file = fileName)

