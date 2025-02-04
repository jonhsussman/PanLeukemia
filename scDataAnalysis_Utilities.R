

## list of functions ####
library(data.table)
library(Matrix)
library(compiler)
library(magrittr)
library(Rtsne)
library(Seurat)
#library(DropletUtils)
#library(cicero)
#library(cisTopic)
#library(ggpubr)
library(RColorBrewer)
#library(heatmaply)
#library(DoubletFinder)
#library(clusterProfiler)
#library(reldist)
#library(MatrixModels)

#library(scABC)
#library(preprocessCore)

#library(chromVAR)
#library(motifmatchr)
#library(SummarizedExperiment)
library(BiocParallel)
#library(JASPAR2016)


# for clustering and/or validation
library(cluster)
#library(factoextra)
#library(fpc)
#library(mclust)

## **** note: for cluster validation, check:  ***** ##
## **** http://www.sthda.com/english/wiki/wiki.php?id_contents=7952  ***##


## do reverse complemente of a DNA sequence
rev.comp <- function(x, rev=TRUE){
  x<-toupper(x)
  y<-rep("N",nchar(x))
  xx<-unlist(strsplit(x,NULL))
  for (bbb in 1:nchar(x))
  {
    if(xx[bbb]=="A") y[bbb]<-"T"    
    if(xx[bbb]=="C") y[bbb]<-"G"    
    if(xx[bbb]=="G") y[bbb]<-"C"    
    if(xx[bbb]=="T") y[bbb]<-"A"
  }
  if(rev==FALSE) 
  {
    for(ccc in (1:nchar(x)))
    {
      if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
    }
  }
  if(rev==T)
  {
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x)))
    {
      zz[ccc]<-y[nchar(x)+1-ccc]
      if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
    }
  }
  return(yy)  
}


segByPeak <- function(peaks){
  peaks[, 'midP' := floor(start/2 + end/2)]
  chrs = unique(peaks$chr)
  domains = NULL
  for(chr0 in chrs){
    bd0 = peaks[chr == chr0]
    setkey(bd0, midP)
    len = nrow(bd0)
    tmp = data.table('chr' = chr0, 'start' = bd0$midP[1:(len - 1)] - 500, 
                     'end' = bd0$midP[2:len] + 500)
    domains = rbind(domains, tmp)
  }
  return(domains)
}
segByPeak = cmpfun(segByPeak)

rebin_matrix2Bin <- function(mtx, resl = 100 * 1000){
  # mtx: matrix wiht rownames as chr-start-end and
  # colnames as cell names
  
  rnames = rownames(mtx)
  mtx_chr = sapply(rnames, function(x) unlist(strsplit(x, '-'))[1])
  chrs = unique(mtx_chr)
  starts = as.numeric(sapply(rnames, function(x) unlist(strsplit(x, '-'))[2]))
  ends = as.numeric(sapply(rnames, function(x) unlist(strsplit(x, '-'))[3]))
  
  rebin_mat = NULL
  for(chr0 in chrs){
    mtx0 = mtx[mtx_chr == chr0, ]
    mtx0 = as.data.table(mtx0)
    mtx0$start = starts[mtx_chr == chr0] 
    mtx0$end = ends[mtx_chr == chr0] 
    
    mtx0[, 'id' := ceiling((start+ end)/resl/2)]
    mtx0[, 'bin_id' := paste0(chr0, '-', id)]
    mtx0[, c('start', 'end', 'id') := NULL]
    rebin_mat = rbind(rebin_mat, mtx0)
  }
  
  rebin_mat = data.table(rebin_mat)
  setkey(rebin_mat, bin_id)
  
  new_mat = rebin_mat[, lapply(.SD, sum), by = bin_id]
  new_mat = new_mat[complete.cases(new_mat)]
  
  feature.names = new_mat$bin_id
  new_mat[, 'bin_id' := NULL]
  
  
  new_mat = as.matrix(new_mat)
  new_mat = as(new_mat, "sparseMatrix")
  rownames(new_mat) = feature.names
  
  return(new_mat)
}
rebin_matrix2Bin = cmpfun(rebin_matrix2Bin)


# filtering of atac matrix
filterMat <- function(atac.mtx, minFrac_in_cell = 0.01, min_depth = 1000,
                      max_depth = 100000){
  depth.cell = Matrix::colSums(atac.mtx)
  atac.mtx = atac.mtx[, depth.cell > min_depth & depth.cell < max_depth]
  frac.in.cell = Matrix::rowSums(atac.mtx > 0)
  atac.mtx = atac.mtx[frac.in.cell > minFrac_in_cell, ]
  return(atac.mtx)
}


# assign gene to nearest peak and mark a gene if its tss within the peak
assignGene2Peak <- function(mtx, gene_ann){
  gene_ann[, 'tss' := ifelse(strand == '+', start, end)]
  peaks = tidyr::separate(data.table(x=rownames(mtx)),
                          col = x,
                          into = c('chr', 'start', 'end'))
  
  
  peaks$peak_name = rownames(mtx)
  class(peaks$start) = 'integer'
  class(peaks$end) = 'integer'
  
  
  
  chrs = unique(peaks$chr)
  peaks_ann = NULL
  for(chr0 in chrs){
    peaks0 = peaks[chr == chr0]
    genes0 = gene_ann[chr == chr0]
    peaks0[, 'id' := which.min(abs(genes0$tss - start/2 - end/2)), by = 'peak_name']
    peaks0[, 'gene_name' := genes0[id, gene_name]]
    peaks0$tss_name = ''
    for(i in 1:nrow(peaks0)){
      tss0 = genes0[tss <= peaks0$end[i] & tss >= peaks0$start[i]]
      if(nrow(tss0) > 0 ) {
        if(peaks0$gene_name[i] %in% tss0$gene_name) peaks0$gene_name[i] <- ''
        peaks0$tss_name[i] = paste(paste0(unique(tss0$gene_name), '-Tss'), 
                                                collapse = ',')
      }
    }
    
    peaks_ann = rbind(peaks_ann, peaks0)
  }
  peaks_ann[, 'id':= NULL]
  
  peaks_ann[, 'peak_new_name' := ifelse(!is.na(gene_name) & nchar(gene_name) > 1, 
                                        paste0(peak_name, ',', gene_name), peak_name)]
  
  
  peaks_ann[, 'peak_new_name' := ifelse(!is.na(tss_name) & nchar(tss_name) > 1, 
                                         paste0(peak_new_name, ',', tss_name), peak_new_name)]
  setkey(peaks_ann, peak_name)
  
  
  
  rownames(mtx) = peaks_ann[rownames(mtx)]$peak_new_name
  
  return(mtx)
  
  
}


# assign gene to nearest peak and mark a gene if its tss within the peak
# input peak_coords with chr-start-end, format
assignGene2Peak_coords <- function(peak_coords, gene_ann){
  gene_ann[, 'tss' := ifelse(strand == '+', start, end)]
  peaks = tidyr::separate(data.table(x = peak_coords),
                          col = x,
                          into = c('chr', 'start', 'end'))
  
  
  peaks$peak_name = peak_coords
  class(peaks$start) = 'integer'
  class(peaks$end) = 'integer'
  
  geneTssInPeak <- function(tss_ids, genes0){
    if(!is.na(tss))
      rr = genes0[tss <= end & tss >= start]$gene_name
    return(paste(rr, collapse = ','))
  }
  
  chrs = unique(peaks$chr)
  peaks_ann = NULL
  for(chr0 in chrs){
    peaks0 = peaks[chr == chr0]
    genes0 = gene_ann[chr == chr0]
    peaks0[, 'id' := which.min(abs(genes0$tss - start/2 - end/2)), by = 'peak_name']
    peaks0[, 'gene_name' := genes0[id, gene_name]]
    peaks0$tss_name = ''
    for(i in 1:nrow(peaks0)){
      tss0 = genes0[tss <= peaks0$end[i] & tss >= peaks0$start[i]]
      if(nrow(tss0) > 0 ) {
        if(peaks0$gene_name[i] %in% tss0$gene_name) peaks0$gene_name[i] <- ''
        peaks0$tss_name[i] = paste(paste0(unique(tss0$gene_name), '-Tss'), 
                                   collapse = ',')
      }
    }
    
    peaks_ann = rbind(peaks_ann, peaks0)
  }
  peaks_ann[, 'id':= NULL]
  
  peaks_ann[, 'peak_new_name' := ifelse(!is.na(gene_name) & nchar(gene_name) > 1, 
                                        paste0(peak_name, ',', gene_name), peak_name)]
  
  
  peaks_ann[, 'peak_new_name' := ifelse(!is.na(tss_name) & nchar(tss_name) > 1, 
                                        paste0(peak_new_name, ',', tss_name), peak_new_name)]
  setkey(peaks_ann, peak_name)
  
  return(peaks_ann[peak_coords, ]$peak_new_name)
  
}


## remove gap between tads
full_seg_tads <- function(tads){
  chrs = unique(tads$chr)
  res = NULL
  for(chr0 in chrs){
    tads0 = tads[chr == chr0]
    bounds = sort(unique(c(tads0$start, tads0$end)))
    len = length(bounds)
    res = rbind(res, data.table('chr' = chr0,
                                'start' = bounds[-len], 'end' = bounds[-1]))
  }
  
  return(res)
}
full_seg_tads = cmpfun(full_seg_tads)


# evaluate clustering given cell*feature matrix object
eval_cluster4mat <- function(mat, cluster.label, alt.cluster.label = NULL, distMethod = 'euclidean'){
 
  Dist = dist(mat, method = distMethod)
  
  label1 = as.numeric(cluster.label)
  if(!is.null(alt.cluster.label)) {
    if(class(alt.cluster.label) == 'character') alt.cluster.label = as.factor(alt.cluster.label)
    alt.cluster.label = as.numeric(alt.cluster.label)
  }
  
  res = cluster.stats(Dist, label1, alt.cluster.label)
  si = silhouette(label1, Dist)
  p <- fviz_silhouette(si, print.summary = F) + 
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  summary.res = list('dunn_single' = res$dunn, 'dunn_ave' = res$dunn2, 
                  'rand' = res$corrected.rand, 'vi' = res$vi,
                  'silhouette' = res$clus.avg.silwidths, 'pl' = p)
  return(summary.res)
}

# evaluate clustering given seurat object
eval_cluster4seurat <- function(seurat.obj, reduction = 'tsne', npc = 30,
                                distMethod = 'euclidean', alt.cluster.label = NULL){
  if(reduction == 'tsne'){
    label1 =  seurat.obj@active.ident
    mat = seurat.obj@reductions$tsne@cell.embeddings
  }
  
  if(reduction == 'umap'){
    label1 =  seurat.obj@active.ident
    mat = seurat.obj@reductions$umap@cell.embeddings
  }
  
  if(reduction == 'pca'){
    label1 =  seurat.obj@active.ident
    mat = seurat.obj@reductions$pca@cell.embeddings[, 1:npc]
  }
  
  Dist = dist(mat, method = distMethod)
  
  label1 = as.numeric(label1)
  if(!is.null(alt.cluster.label)) {
    if(class(alt.cluster.label) == 'character') alt.cluster.label = as.factor(alt.cluster.label)
    alt.cluster.label = as.numeric(alt.cluster.label)
  }
  
  res = cluster.stats(Dist, label1, alt.cluster.label)
  si = silhouette(label1, Dist)
  p <- fviz_silhouette(si, print.summary = F) + 
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  summary.res = list('dunn_single' = res$dunn, 'dunn_ave' = res$dunn2, 
                  'rand' = res$corrected.rand, 'vi' = res$vi,
                  'silhouette' = res$clus.avg.silwidths, 'pl' = p)
  return(summary.res)
}

# do normalization, pca using Seurat
doBasicSeurat <- function(mtx, npc = 50, top.variable = 0.2, doLog = T, 
                           doScale = T, doCenter = T, assay = 'ATAC', 
                           reg.var = 'nCount_ATAC'){
  
 # top.variabl -- use top most variable features
  if(doLog) mtx = round(log1p(mtx) / log(2))
  seurat.obj = CreateSeuratObject(mtx, project = 'scATAC', assay = assay,
                                  names.delim = '-')
  cell.names = colnames(mtx)
  
  #seurat.obj <- NormalizeData(seurat.obj, normalization.method = 'LogNormalize',
  #                            scale.factor = 1e4)
  
  seurat.obj <- FindVariableFeatures(object = seurat.obj, 
                                     selection.method = 'vst', 
                                     nfeatures = floor(nrow(mtx) * top.variable))
  seurat.obj <- ScaleData(object = seurat.obj, 
                          features = VariableFeatures(seurat.obj), 
                          vars.to.regress = reg.var, do.scale = doScale, do.center = doCenter)
  
  
  seurat.obj <- RunPCA(object = seurat.obj, 
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  
  return(seurat.obj)
}
doBasicSeurat = cmpfun(doBasicSeurat)


regress_on_pca <- function(seurat.obj, reg.var = 'nCount_ATAC'){
  
  pcs = seurat.obj@reductions$pca@cell.embeddings
  pcs.reg = pcs
  for(i in 1:length(reg.var)){
    
    reg.var0 = seurat.obj[[reg.var[i]]][[1]]
    pcs.reg = apply(pcs.reg, 2, function(x) lm(x ~ reg.var0)$residual )
    
  }
   colnames(pcs.reg) = colnames(pcs)
  seurat.obj@reductions$pca@cell.embeddings = pcs.reg
  return(seurat.obj)
}

# do normalization using log, tf-idf, or none, regress out confounds on pca or not 
doBasicSeurat_atac <- function(mtx, npc = 50, top.variable = 0.2, 
                               norm_by = c('log', 'tf-idf', 'none'),
                              doScale = T, doCenter = T, assay = 'ATAC',
                              reg.var = 'nCount_ATAC', regressOnPca = T,
                              project = 'scATAC'){
  
  # top.variabl -- use top most variable features
  seurat.obj = CreateSeuratObject(mtx, project = project, assay = assay,
                                  names.delim = '-')
 
  if(norm_by == 'log') seurat.obj@assays[[assay]]@data <- log1p(mtx) / log(2)
  if(norm_by == 'tf-idf') seurat.obj@assays[[assay]]@data <- TF.IDF(mtx, verbose = F)
  
  seurat.obj <- FindVariableFeatures(object = seurat.obj,
                                     selection.method = 'vst',
                                     nfeatures = floor(nrow(mtx) * top.variable))
  
  if(regressOnPca){
    reg.var0 = NULL
  }else{
    reg.var0 = reg.var
  }
  
  seurat.obj <- ScaleData(object = seurat.obj,
                          features = VariableFeatures(seurat.obj),
                          vars.to.regress = reg.var0, do.scale = doScale,
                          do.center = doCenter)
  
  
  #seurat.obj <- RunPCA(object = seurat.obj,
  #                     features = VariableFeatures(object = seurat.obj),
  #                     verbose = FALSE, seed.use = 10, npc = npc)
  seurat.obj <- RunPCA(object = seurat.obj,
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, npc = npc)
  if(length(reg.var) > 0 & regressOnPca) seurat.obj = regress_on_pca(seurat.obj, reg.var)
  
  return(seurat.obj)
}

doBasicSeurat_atac = cmpfun(doBasicSeurat_atac)

doBasicSeurat_atac_updated <- function(mtx, npc = 30, top.variable = 5000, 
                                       norm_by = c('log', 'tf-idf', 'none'),
                                       doScale = T, doCenter = T, assay = 'ATAC',
                                       reg.var = 'nCount_ATAC', regressOnPca = T,
                                       project = 'scATAC', vap.min.frac = 0,
                                       meta.data = NULL){
  
  # top.variabl -- use top most variable features
  seurat.obj = CreateSeuratObject(mtx, project = project, assay = assay,
                                  names.delim = '-',
                                  meta.data = meta.data)
  
  if(norm_by == 'log') seurat.obj@assays[[assay]]@data <- log1p(mtx) / log(2)
  if(norm_by == 'tf-idf') seurat.obj@assays[[assay]]@data <- TF.IDF(mtx, verbose = F)
  
  nvap = ifelse(top.variable > 1, top.variable, floor(top.variable * ncol(mtx)))
  seurat.obj <- FindVariableFeatures(object = seurat.obj,
                                     selection.method = 'vst',
                                     nfeatures = nvap)
  vaps = VariableFeatures(seurat.obj)
  peak.frac = rowMeans(mtx > 0)
  excludePks.fromVAP = names(which(peak.frac < vap.min.frac))
  vaps = setdiff(vaps, excludePks.fromVAP)
  
  if(length(vaps) < 10) stop('Top few VAPs left!')
  ## redo normalization using vap
  if(norm_by == 'tf-idf'){
    mtx.norm = TF.IDF(mtx[vaps, ])
    tmp <- mtx[setdiff(rownames(mtx), vaps), ]
    data0 <- rbind(mtx.norm, tmp)
    seurat.obj[[assay]]@data = data0[rownames(mtx), ]
    rm(data0, tmp, mtx.norm)
  }
  
  
  if(regressOnPca){
    reg.var0 = NULL
  }else{
    reg.var0 = reg.var
  }
  VariableFeatures(seurat.obj) <- vaps
  seurat.obj <- ScaleData(object = seurat.obj,
                          features = vaps,
                          vars.to.regress = reg.var0, do.scale = doScale,
                          do.center = doCenter)
  
  seurat.obj <- RunPCA(object = seurat.obj,
                       features = vaps,
                       verbose = FALSE, npc = npc)
  if(length(reg.var) > 0 & regressOnPca) seurat.obj = regress_on_pca(seurat.obj, reg.var)
  
  return(seurat.obj)
}

doBasicSeurat_atac_updated = cmpfun(doBasicSeurat_atac_updated)

runSeurat_Atac <- function(mtx, npc = 50, top_variable_features = 0.2, 
                           doScale = T, doCenter = T, assay = 'ATAC',
                           reg.var = NULL, norm_by = 'log', project = 'scATAC'){
  
  # top.variabl -- use top most variable features
  seurat.obj = CreateSeuratObject(mtx, project = project, assay = assay,
                                  names.delim = '-', min.cells = 1,
                                  min.features = 1)
  if(norm_by == 'log') seurat.obj[[assay]]@data <- log1p(seurat.obj[[assay]]@counts)/log(2)
  if(norm_by == 'tf-idf') seurat.obj[[assay]]@data <- TF.IDF(seurat.obj[[assay]]@counts)
  nvap = ifelse(top_variable_features > 1, top_variable_features, 
                floor(nrow(mtx) * top_variable_features))
  
  seurat.obj <- FindVariableFeatures(object = seurat.obj,
                                     selection.method = 'vst',
                                     nfeatures = nvap)
  
  ## remove variable features only accessible in less than 1% of cells
  mtx = seurat.obj[[assay]]@counts
  rs = Matrix::rowMeans(mtx > 0)
  rare.features = names(which(rs < 0.01))
  vaps = VariableFeatures(seurat.obj)
  vaps = setdiff(vaps, rare.features)
  niter = 0
  while(length(vaps) < 500 & nvap > 500){
    niter = niter + 1
    nvap = nvap + 2000
    seurat.obj <- FindVariableFeatures(object = seurat.obj,
                                       selection.method = 'vst',
                                       nfeatures = min(nvap, nrow(seurat.obj)))
    vaps = VariableFeatures(seurat.obj)
    vaps = setdiff(vaps, rare.features)
    if(niter > 5) stop('Too many variable features were filtered, 
                       please specify a large Top_Variable_Features
                       in the configure file!')
  }
  VariableFeatures(seurat.obj) <- vaps
  
  ## redo normalization using vap if norm by tf-idf
  if(norm_by == 'tf-idf'){
    mtx.norm = TF.IDF(mtx[vaps, ])
    tmp <- mtx[setdiff(rownames(mtx), vaps), ]
    data0 <- rbind(mtx.norm, tmp)
    seurat.obj[[assay]]@data = data0[rownames(mtx), ]
    rm(data0, tmp, mtx.norm)
  }
  
  seurat.obj <- ScaleData(object = seurat.obj,
                          features = VariableFeatures(seurat.obj),
                          vars.to.regress = NULL, do.scale = doScale,
                          do.center = doCenter)
  
  
  seurat.obj <- RunPCA(object = seurat.obj,
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  if(length(reg.var) > 0 ) seurat.obj = regress_on_pca(seurat.obj, reg.var)
  
  
  return(seurat.obj)
  }
runSeurat_Atac = cmpfun(runSeurat_Atac)

# do normalization, pca using Seurat
doBasicSeurat_RNA <- function(mtx, npc = 50, top.variable = 0.2, pmito.upper = 0.2,
                           doScale = T, doCenter = T, min_nCount = 1500,
                           max_nCount = 15000, reg.var = 'nCount_RNA', sct=F){

 ## top.variabl -- use top most variable features

 # filter cells with high percentage of mitocondria genes

  nCount = Matrix::colSums(mtx)
  mtx = mtx[, nCount < max_nCount & nCount > min_nCount]
  
  mito.features <- grep(pattern = "^MT-", 
                      x = rownames(x = mtx), value = TRUE)

  perc.mito = Matrix::colSums(mtx[mito.features, ])/Matrix::colSums(mtx)

  mtx = mtx[, perc.mito <= pmito.upper]
  perc.mito = perc.mito[perc.mito <= pmito.upper]


 # create seurat object
  seurat.obj = CreateSeuratObject(mtx, project = 'scRNA', assay = 'RNA',
                                  names.delim = '-', min.cells = 0, min.features = 0)
  

  # add perc.mito to seurat objects
  cnames = colnames(seurat.obj)
  tmp.mito = data.table('perc' = perc.mito, 'cname' = names(perc.mito))
  setkey(tmp.mito, cname)
  seurat.obj@meta.data[['perc.mito']] = tmp.mito[cnames]$perc
  
  
  seurat.obj <- subset(x = seurat.obj, subset = (nFeature_RNA < 10000))
  
  if(!sct){
    seurat.obj <- NormalizeData(seurat.obj, normalization.method = 'LogNormalize',
                                scale.factor = 1e4)
    seurat.obj <- FindVariableFeatures(object = seurat.obj, 
                                       selection.method = 'vst', 
                                       nfeatures = floor(nrow(mtx) * top.variable))
    seurat.obj <- ScaleData(object = seurat.obj, 
                            features = VariableFeatures(seurat.obj), 
                            vars.to.regress = reg.var, do.scale = doScale, do.center = doCenter)
    
    
  }else{
    seurat.obj <- SCTransform(seurat.obj, vars.to.regress = reg.var, verbose = F,
                              variable.features.n = floor(nrow(mtx) * top.variable))
  }
  
  #seurat.obj <- RunPCA(object = seurat.obj, 
  #                     features = VariableFeatures(object = seurat.obj),
  #                     verbose = FALSE, seed.use = 10, npc = npc)
  seurat.obj <- RunPCA(object = seurat.obj, 
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE,  npc = npc)
  
  return(seurat.obj)
}
doBasicSeurat_RNA = cmpfun(doBasicSeurat_RNA)



# Find doublets
FindDoublets <- function(seurat.rna, PCs = 1:50, exp_rate = 0.02, sct = FALSE){
  # sct--do SCTransform or not
  
  ## pK identification
  sweep.res.list <- paramSweep_v3(seurat.rna, PCs = PCs, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  ## Homotypic Doublet proportion Estimate
  annotations <- seurat.rna@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(exp_rate * length(seurat.rna$seurat_clusters))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  seurat.rna <- doubletFinder_v3(seurat.rna, PCs = PCs, pN = 0.25,
                                 pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, 
                                 sct = sct)
  
  seurat.rna <- doubletFinder_v3(seurat.rna, PCs = PCs, pN = 0.25, 
                                 pK = 0.09, nExp = nExp_poi.adj,
                                 reuse.pANN = paste0("pANN_0.25_0.09_", nExp_poi), 
                                 sct = sct)
  doublet_var = paste0('DF.classifications_0.25_0.09_', nExp_poi.adj)
  seurat.rna[['Doublet_Singlet']] = seurat.rna[[doublet_var]]
  
  mnames = names(seurat.rna@meta.data)
  seurat.rna@meta.data[, grep(mnames, pattern = '0.25_0.09')] <- NULL
  #seurat.rna = subset(seurat.rna, Doublet_Singlet == 'Singlet')
  return(seurat.rna)
}


doSeurat_rmDoublets <- function(sampleID, exp_rate = 0.05, pmito.upper = 0.15,
                                min_nCount = 1500, max_nCount = 50000,
                                clusterOn = 'pca', npc = 50,
                                resolution = 0.5){
  
  dir0 = '/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/CellRangerResults/June19_2019/'
  
  sampleName = paste0('MLL_', sampleID, '_scRNA')
  mtx = Read10X(paste0(dir0, sampleName, '/', sampleName, '/outs/filtered_feature_bc_matrix'))
  colnames(mtx) = paste0(sampleID, '_', colnames(mtx))
  
  # remove red blood cells
  if('HBB' %in% rownames(mtx)){
    hbb_exp = mtx['HBB', ]
    mtx = mtx[, hbb_exp < 3]
  }
  
  
  ##try remove doublets instead of manully filter cells with larger UMI
  seurat.obj = doBasicSeurat_RNA(mtx, min_nCount = min_nCount, max_nCount = max_nCount,
                                 npc = npc, pmito.upper = pmito.upper)
  seurat.obj = FindNeighbors(seurat.obj, reduction = clusterOn, dims = 1:npc)
  seurat.obj = FindClusters(seurat.obj, resolution = resolution)
  seurat.obj = FindDoublets(seurat.obj, exp_rate = exp_rate)
  seurat.obj = RunUMAP(seurat.obj, dims = 1:npc)
  p1 <- DimPlot(seurat.obj, reduction = 'umap', 
                group.by = 'Doublet_Singlet') + ggtitle('With Doublets')
  
  ## remove doublets and do PCA and clustering again
  seurat.obj = subset(seurat.obj, Doublet_Singlet == 'Singlet')
  seurat.obj = doBasicSeurat_RNA(seurat.obj@assays$RNA@counts, 
                                 min_nCount = min_nCount, max_nCount = max_nCount)
  seurat.obj = RunUMAP(seurat.obj, dims = 1:npc)
  seurat.obj = RunTSNE(seurat.obj, dims = 1:npc)
  seurat.obj = FindNeighbors(seurat.obj, reduction = 'pca', dims = 1:npc)
  seurat.obj = FindClusters(seurat.obj, resolution = resolution)
  
  p2 <- DimPlot(seurat.obj, reduction = 'umap', label = T) + ggtitle('Doublets Removed')
  
  
  
  
  ## load to viscello
  inputs = prepInput4Cello(mtx = seurat.obj@assays$RNA@counts, 
                           seurat.obj = seurat.obj,
                           cello.name = sampleName,
                           assay = 'RNA', 
                           extraDR = T, cluster_excluded = NULL)
  
  # save inputs for future use
  celloInputDir = paste0('Input4VisCello/', sampleName)
  system(paste0('mkdir -p ', celloInputDir))
  saveRDS(inputs$eset, file = paste0(celloInputDir, '/eset.rds'))
  saveRDS(inputs$clist, file = paste0(celloInputDir, '/clist.rds'))
  
  ggsave(CombinePlots(plots = list(p1, p2)), device = 'eps', width = 14, 
         height = 6, filename = paste0('Figures/scRNA/MLL_', 
                                       sampleID, '/', sampleName, '_umap_with_doublets.eps'))
  
  
  saveRDS(seurat.obj, paste0('Seurat_Objects/scRNA/seurat_', sampleName, '_doubletRemoved.rds'))
  
  return(seurat.obj)
}

##bcPrefix can be set as sampleName or ID
doSeurat_rmDoublets_dir <- function(mtx, exp_rate = 0.05, pmito.upper = 0.15,
                                min_nCount = 1500, max_nCount = 50000,
                                clusterOn = 'pca', npc = 50,
                                resolution = 0.5, bcPrefix = 'pbmc'){
  
  
  #mtx = Read10X(filtered_mtx_dir)
  colnames(mtx) = paste0(bcPrefix, '_', colnames(mtx))
  
  # remove red blood cells
  if('HBB' %in% rownames(mtx)){
    hbb_exp = mtx['HBB', ]
    mtx = mtx[, hbb_exp < 3]
  }
  
  
  ##try remove doublets instead of manully filter cells with larger UMI
  seurat.obj = doBasicSeurat_RNA(mtx, min_nCount = min_nCount, max_nCount = max_nCount,
                                 npc = npc, pmito.upper = pmito.upper)
  #seurat.obj = FindNeighbors(seurat.obj, reduction = clusterOn, dims = 1:npc)
  #seurat.obj = FindClusters(seurat.obj, resolution = resolution)
  #seurat.obj = FindDoublets(seurat.obj, exp_rate = exp_rate)
  #seurat.obj = RunUMAP(seurat.obj, dims = 1:npc)
  #p1 <- DimPlot(seurat.obj, reduction = 'umap', 
  #              group.by = 'Doublet_Singlet') + ggtitle('With Doublets')
  
  ## remove doublets and do PCA and clustering again
  #seurat.obj = subset(seurat.obj, Doublet_Singlet == 'Singlet')
  #seurat.obj = doBasicSeurat_RNA(seurat.obj@assays$RNA@counts, 
  #                               min_nCount = min_nCount, max_nCount = max_nCount)
  seurat.obj = RunUMAP(seurat.obj, dims = 1:npc)
  seurat.obj = RunTSNE(seurat.obj, dims = 1:npc)
  seurat.obj = FindNeighbors(seurat.obj, reduction = 'pca', dims = 1:npc)
  seurat.obj = FindClusters(seurat.obj, resolution = resolution)
  #p2 <- DimPlot(seurat.obj, reduction = 'umap', label = T) + 
  #  ggtitle('Doublets Removed')
  
  #ggsave(CombinePlots(plots = list(p1, p2)), device = 'eps', width = 14, 
  #       height = 6, filename = paste0('Figures/scRNA/', bcPrefix, '_umap_with_doublets.eps'))
  
  seurat.obj$sample = bcPrefix
  return(seurat.obj)
}

doSeurat_rmDoublets_dir_new <- function(filtered_mtx_dir, exp_rate = 0.05, pmito.upper = 0.15,
                                    min_nCount = 1500, max_nCount = 50000,
                                    clusterOn = 'pca', npc = 50,
                                    resolution = 0.5, bcPrefix = 'pbmc'){
  
  
  mtx = Read10X(filtered_mtx_dir)
  colnames(mtx) = paste0(bcPrefix, '_', colnames(mtx))
  
  # remove red blood cells
  if('HBB' %in% rownames(mtx)){
    hbb_exp = mtx['HBB', ]
    mtx = mtx[, hbb_exp < 3]
  }
  
  
  ##try remove doublets instead of manully filter cells with larger UMI
  seurat.obj = doBasicSeurat_RNA(mtx, min_nCount = min_nCount, max_nCount = max_nCount,
                                 npc = npc, pmito.upper = pmito.upper)
  #seurat.obj = FindNeighbors(seurat.obj, reduction = clusterOn, dims = 1:npc)
  #seurat.obj = FindClusters(seurat.obj, resolution = resolution)
  #seurat.obj = FindDoublets(seurat.obj, exp_rate = exp_rate)
  #seurat.obj = RunUMAP(seurat.obj, dims = 1:npc)
  #p1 <- DimPlot(seurat.obj, reduction = 'umap', 
  #              group.by = 'Doublet_Singlet') + ggtitle('With Doublets')
  
  ## remove doublets and do PCA and clustering again
  #seurat.obj = subset(seurat.obj, Doublet_Singlet == 'Singlet')
  #seurat.obj = doBasicSeurat_RNA(seurat.obj@assays$RNA@counts, 
  #                               min_nCount = min_nCount, max_nCount = max_nCount)
  seurat.obj = RunUMAP(seurat.obj, dims = 1:npc)
  seurat.obj = RunTSNE(seurat.obj, dims = 1:npc)
  seurat.obj = FindNeighbors(seurat.obj, reduction = 'pca', dims = 1:npc)
  seurat.obj = FindClusters(seurat.obj, resolution = resolution)
  #p2 <- DimPlot(seurat.obj, reduction = 'umap', label = T) + 
  #  ggtitle('Doublets Removed')
  
  #ggsave(CombinePlots(plots = list(p1, p2)), device = 'eps', width = 14, 
  #       height = 6, filename = paste0('Figures/scRNA/', bcPrefix, '_umap_with_doublets.eps'))
  
  seurat.obj$sample = bcPrefix
  return(seurat.obj)
}


# basic seurat plot
basicSeuratPlot <- function(org.seurat.obj){
  p1 <- VizDimLoadings(object = org.seurat.obj, dims = 1:2, nfeatures = 20)
  
  p2 <- DimHeatmap(object = org.seurat.obj, dims = 1:6, cells = 100,   balanced = TRUE)
  p3 <- ElbowPlot(object = org.seurat.obj, ndims = 50)
  mfeatures = names(org.seurat.obj@meta.data)
  pfeatures = mfeatures[grepl(mfeatures, pattern = 'nCount_|mito')]
  p4 <- VlnPlot(org.seurat.obj, pfeatures)
  return(list(p1, p2, p3, p4))
}


# integrate analysis
# data were given in unit of TAD
doSeurat_integrate_tad <- function(atac.mtx, gene.mtx, qLorm = T){
  colnames(atac.mtx) = paste0('atac-', colnames(atac.mtx))
  comb.mtx = cbind(gene.mtx, atac.mtx)
  dtype = sapply(colnames(comb.mtx), function(x) ifelse(grepl(x, pattern = '^atac-'), 'ATAC', 'RNA'))
  
  if(qLorm) {
    rnames = rownames(comb.mtx)
    cnames = colnames(comb.mtx)
    rand.mat = matrix(runif(nrow(comb.mtx) * ncol(comb.mtx), 0, 10^(-6)), nrow(comb.mtx))
    comb.mtx <- normalize.quantiles(as.matrix(comb.mtx) + rand.mat)
    colnames(comb.mtx) = cnames
    rownames(comb.mtx) = rnames
  }
  seurat.integrate <- CreateSeuratObject(comb.mtx, assay = 'Comb')
  seurat.integrate@meta.data$orig.ident = dtype

  #seurat.obj <- NormalizeData(seurat.integrate, normalization.method = 'LogNormalize',
  #                            scale.factor = 1e4)
  
  seurat.integrate <- FindVariableFeatures(object = seurat.integrate, 
                                     selection.method = 'dispersion', 
                                     nfeatures = nrow(comb.mtx))
  seurat.integrate <- ScaleData(object = seurat.integrate, 
                          features = VariableFeatures(seurat.integrate), 
                          vars.to.regress = c('orig.ident', 'nCount_Comb'), do.scale = F, do.center = F)
  
  
  seurat.integrate <- RunPCA(object = seurat.integrate, 
                       features = VariableFeatures(object = seurat.integrate),
                       verbose = FALSE, seed.use = 10, npc = 50)
  return(seurat.integrate)
}


## map gene to atac peak
gene2peak <- function(gene_set, peaks, gene_ann){
  # should include tss information in gene_list
  gene_list = gene_ann[gene_name %in% gene_set, ]
  chrs = unique(gene_list$chr)
  gene_new = NULL
  peaks[, 'midP' := start/2 + end/2]
  for(chr0 in chrs){
    gene0 = gene_list[chr == chr0, ]
    peaks0 = peaks[chr == chr0]
    gene0[, 'peak_id' := which(tss >= peaks0$start & tss <= peaks0$end), by = gene_id]
    gene0[, 'peak_id' := ifelse(is.na(peak_id), which.min(abs(tss - peaks0$midP)), peak_id), by = gene_name]
    gene0[, 'peak' := peaks0[peak_id]$pos]
    
    gene_new = rbind(gene_new, gene0)
  }
  return(gene_new)
}


# plot atac signal around Tss, given a gene set and seurat atac object
plotSignal_Tss <- function(genes, org.seurat.obj, reduction = 'tsne'){
  ## genes should include gene_name and the corresponding peak information
  tmp_plot = list()
  for(i in 1:nrow(genes)){
    tmp_plot[[i]] <- FeaturePlot(object = org.seurat.obj, reduction = reduction,  feature = genes$peak[i]) + labs(title = genes$gene_name[i])
  } 
  
  pp <- ggarrange(plotlist = tmp_plot, nrow = 2, ncol = ceiling(length(tmp_plot)/2))
  return(pp)
}

read10X_ATAC <- function(dirt){
  mtx_path <- paste0(dirt, "matrix.mtx")
  feature_path <- paste0(dirt, "peaks.bed")
  barcode_path <- paste0(dirt, "barcodes.tsv")
  
  
  features <-readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature, sep = '-')
  barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)
  
  mtx <-  Matrix::readMM(mtx_path) %>%
    magrittr::set_rownames(features$feature)%>%
    magrittr::set_colnames(barcodes$barcode) 
  
  return(mtx)
}

read_mtx_scATACpro <- function(mtx_path){
  #mtx_path <- paste0(dirt, "matrix.mtx")
  mtx.dir = dirname(mtx_path)
  feature_path <- paste0(mtx.dir, "/features.txt")
  barcode_path <- paste0(mtx.dir, "/barcodes.txt")
  
  
  features <- fread(feature_path, header = F)
  barcodes <- fread(barcode_path, header = F)
  
  mtx <-  Matrix::readMM(mtx_path) %>%
    magrittr::set_rownames(features$V1)%>%
    magrittr::set_colnames(barcodes$V1)
  
  return(mtx)
}


compare_cellCall4ATAC <- function(cellRanger_dir, fdr = 0.001, lower = 500, upper = NULL){
  filter_dir = paste0(cellRanger_dir, 'filtered_peak_bc_matrix/')
  raw_dir = paste0(cellRanger_dir, 'raw_peak_bc_matrix/')
  
  mat <- read10X_ATAC(raw_dir)
  
  
  set.seed(2019)
  cell.out <- emptyDrops(mat, lower = lower, retain = upper)
  
  filter.out <- cell.out[complete.cases(cell.out), ]
  
  #is.cell <- (cell.out$FDR <= fdr)
  #plot(cell.out$Total, -cell.out$LogProb, col=ifelse(is.cell, "red", "black"),
  #     xlab="Total count", ylab="-Log Probability", log = 'x')
  
  
  filter.out = filter.out[filter.out$FDR <= fdr, ]
  
  rm(mat)
  
  
  mat <- read10X_ATAC(filter_dir)
  
  overlapped.cell <- intersect(rownames(filter.out), colnames(mat))
  dim(mat)
  perc.in.cellranger <- length(overlapped.cell)/nrow(filter.out)
  perc.in.emptydrop <- length(overlapped.cell)/ncol(mat)
  
  
  
  output = list('emptyDrop_cells' = filter.out, 'perc.in.cellranger' = round(perc.in.cellranger, 3),
                'perc.in.emptydrop' = round(perc.in.emptydrop, 3))
  
  return(output)
}


# do cicero given a Seurat object
doCicero_gascore <- function(seurat.obj, reduction = 'umap', chr_sizes,
                     gene_ann, coaccess_thr = 0.25, peak_in_cell = 0.025){
  ## gene_ann: the first four columns: chr, start, end, gene name
  set.seed(2019) 

  mtx = GetAssayData(seurat.obj, slot = 'counts')
  # change rownames using _ to delimited
  rnames = rownames(mtx)
  new.rnames = sapply(rnames, function(x) unlist(strsplit(x, ','))[1])
  new.rnames = sapply(new.rnames, function(x) gsub('-', '_', x))
  rownames(mtx) <- new.rnames
  
  #filter some peaks
  if(peak_in_cell > 0){
    rr = Matrix::rowMeans(mtx > 0)
    mtx = mtx[rr > peak_in_cell, ]
  }
  
  dt = reshape2::melt(as.matrix(mtx), value.name = 'count')
  rm(mtx)
  dt = dt[dt$count > 0, ]
  input_cds <- make_atac_cds(dt, binarize = T)
  rm(dt)
  input_cds <- detectGenes(input_cds)
  input_cds <- estimateSizeFactors(input_cds)
  
  if(reduction == 'tsne') {
    redu.coords = seurat.obj@reductions$tsne@cell.embeddings
  }
  if(reduction == 'umap') {
    redu.coords = seurat.obj@reductions$umap@cell.embeddings
  }
  
  #make the cell id consistet
  
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = redu.coords)
  
  ## get connections
  
  conns <- run_cicero(cicero_cds, chr_sizes) 
  
  ## get cicero gene activity score
  names(gene_ann)[4] <- "gene"
  
  input_cds <- annotate_cds_by_site(input_cds, gene_ann)
  
  # generate unnormalized gene activity matrix
  unnorm_ga <- build_gene_activity_matrix(input_cds, conns)
  
  # make a list of num_genes_expressed
  num_genes <- pData(input_cds)$num_genes_expressed
  names(num_genes) <- row.names(pData(input_cds))
  
  # normalize
  cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
  
  # if you had two datasets to normalize, you would pass both:
  # num_genes should then include all cells from both sets
  #unnorm_ga2 <- unnorm_ga
  #cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga, unnorm_ga2), num_genes)
  
  conns = data.table(conns)
  conns = conns[coaccess > coaccess_thr, ]
  res = list('conns' = conns, 'ga_score' = cicero_gene_activities)
  return(res)
}


# do cicero given a Seurat object, just return the connection 
doCicero_conn <- function(seurat.obj, reduction = 'tsne', chr_sizes, npc = 30){
  ## gene_ann: the first four columns: chr, start, end, gene name
  set.seed(2019) 

  mtx = GetAssayData(seurat.obj, slot = 'counts')
  # change rownames using _ to delimited
  rnames = rownames(mtx)
  new.rnames = sapply(rnames, function(x) gsub('-', '_', x))
  rownames(mtx) <- new.rnames
  
  dt = reshape2::melt(as.matrix(mtx), value.name = 'count')
  rm(mtx)
  dt = dt[dt$count > 0, ]
  input_cds <- make_atac_cds(dt, binarize = T)
  rm(dt)
  input_cds <- detectGenes(input_cds)
  input_cds <- estimateSizeFactors(input_cds)
  
  if(reduction == 'tsne') {
    if(is.null(seurat.obj@reductions$tsne))
      seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc)
    redu.coords = seurat.obj@reductions$tsne@cell.embeddings
  }
  if(reduction == 'umap') {
    if(is.null(seurat.obj@reductions$umap))
      seurat.object <- RunUMAP(seurat.object, dims = 1:npc)
    redu.coords = seurat.obj@reductions$umap@cell.embeddings
  }
  
  #make the cell id consistet
  
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = redu.coords)
  
  ## get connections
  
  conns <- run_cicero(cicero_cds, chr_sizes) 
  

  return(conns)
}



## query the resoltuion parameters given a seurat object and the number of clusters
## using binary seach
queryResolution4Seurat <- function(seurat.obj, k = 10, reduction = 'umap', npc = 20, 
                      min_resl = 0.1, max_resl = 1, max_iter = 15, doPCA = F){
  max.dim = ifelse(reduction == 'pca', npc, 2)
  if(doPCA) {
    seurat.obj <- RunPCA(seurat.obj, npcs = npc, verbose = F)
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc, verbose = F)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:npc, verbose = F)
  }

  
  seurat.obj <- FindNeighbors(seurat.obj, reduction = reduction, verbose = F, dims = 1:max.dim)
  tmp.cluster1 <- FindClusters(seurat.obj, resolution = min_resl)@active.ident
  tmp.cluster2 <- FindClusters(seurat.obj, resolution = max_resl)@active.ident
  
  
  

  len1 = length(levels(tmp.cluster1))
  len2 = length(levels(tmp.cluster2))

  k1 = k2 = 0
  while(len1 > k ){
   
    k1 = k1 + 1
    message('min_resl too large, trying to divided it by  2')
    min_resl = min_resl/2
    tmp.cluster1 <- FindClusters(seurat.obj, resolution = min_resl)@active.ident
    len1 = length(levels(tmp.cluster1))
    if(k1 == 10) stop('Please specify a much smaller min_res')
  }

  while(len2 < k){
    k2 = k2 + 1
    message('max_resl too small, trying to multiply it by 2')
    max_resl = max_resl * 2
    tmp.cluster2 <- FindClusters(seurat.obj, resolution = max_resl)@active.ident
    len2 = length(levels(tmp.cluster2))
    if(k2 == 10) stop('Please specify a much bigger max_res')
  }
  if(len1 == k) {
    return(min_resl)
  }

  if(len2 == k) {
    return(max_resl)
  }

  # repeat in other case
  
  i = 0
  repeat{
    i = i + 1
    resl0 = min_resl/2 + max_resl/2
    
    tmp.cluster <- FindClusters(seurat.obj, resolution = resl0)@active.ident
      
    len = length(levels(tmp.cluster)) 
    if(len == k){
      return(resl0)
    }
    if(len < k){
      min_resl = resl0
      len1 = len
    }
    if(len > k){
      max_resl = resl0
      len2 = len
    }
    if(i == max_iter) break
  }
  return(resl0)
}
queryResolution4Seurat = cmpfun(queryResolution4Seurat)



## query the resoltuion parameters given a seurat object (transformed from cistopic object) and the number of clusters
## using binary seach
queryResolution4Topic <- function(seurat.obj, k = 10,  min_resl = 0.1, max_resl = 1, max_iter = 15){

  
  # skip find neighbors 
 
  tmp.cluster1 <- FindClusters(seurat.obj, resolution = min_resl, graph.name = 'snn')@active.ident
  tmp.cluster2 <- FindClusters(seurat.obj, resolution = max_resl, graph.name = 'snn')@active.ident
  
  
  

  len1 = length(levels(tmp.cluster1))
  len2 = length(levels(tmp.cluster2))

  k1 = k2 = 0
  while(len1 > k ){
   
    k1 = k1 + 1
    message('min_resl too large, trying to divided it by  2')
    min_resl = min_resl/2
    tmp.cluster1 <- FindClusters(seurat.obj, resolution = min_resl, graph.name = 'snn')@active.ident
    len1 = length(levels(tmp.cluster1))
    if(k1 == 10) stop('Please specify a much smaller min_res')
  }

  while(len2 < k){
    k2 = k2 + 1
    message('max_resl too small, trying to multiply it by 2')
    max_resl = max_resl * 2
    tmp.cluster2 <- FindClusters(seurat.obj, resolution = max_resl, graph.name = 'snn')@active.ident
    len2 = length(levels(tmp.cluster2))
    if(k2 == 10) stop('Please specify a much bigger max_res')
  }


  if(len1 == k) {
    return(min_resl)
  }

  if(len2 == k) {
    return(max_resl)
  }

  # repeat in other case
  
  i = 0
  repeat{
    i = i + 1
    resl0 = min_resl/2 + max_resl/2
    
    tmp.cluster <- FindClusters(seurat.obj, resolution = resl0, graph.name = 'snn')@active.ident
      
    len = length(levels(tmp.cluster)) 
    if(len == k){
      return(resl0)
    }
    if(len < k){
      min_resl = resl0
      len1 = len
    }
    if(len > k){
      max_resl = resl0
      len2 = len
    }
    if(i == max_iter) break
  }
  return(resl0)
}
queryResolution4Topic = cmpfun(queryResolution4Topic)



# plot cicero gene activity score, given a gene set, cicero gene activity score and seurat atac object
plotCicero_ascore <- function(genes, cicero_ascore, seurat.atac, reduction = 'tsne'){
  ## genes should include gene_name and gene_id information 
  ## gene_id is the rownames of cicero_ascore, cell id is the colnames of cicero_ascore
  
  cicero_ascore = cicero_ascore[rownames(cicero_ascore) %in% genes$gene_id, ]
  
  if(nrow(cicero_ascore) == 0) stop('No activity score found for these genes!')
  genes = genes[gene_id %in% rownames(cicero_ascore), ]

  # add the gene activity score as a metadata feature
  if(any(colnames(seurat.atac) != colnames(cicero_ascore))) stop('Cell id not consistent!')

  for(gene_id0 in genes$gene_id){
    seurat.atac@meta.data[[gene_id0]] <- cicero_ascore[gene_id0, ]
  }

  tmp_plot = list()
  for(i in 1:nrow(genes)){
    tmp_plot[[i]] <- FeaturePlot(object = seurat.atac, reduction = reduction,  feature = genes$gene_id[i]) + 
    labs(title = genes$gene_name[i])
  } 
  
  pp <- ggarrange(plotlist = tmp_plot, ncol = 2, nrow = ceiling(length(tmp_plot)/2))
  return(pp)
}


basicCluster <- function(reduced.mtx, method = 'kmeans', k = 5){
  if(method == 'kmeans'){
    res = kmeans(reduced.mtx, centers = k)
    cl.label = res$cluster
  }

  if(method == 'hclust'){
    if(is.null(k)) stop('Need specify k: the number of cluster')
    d <- dist(reduced.mtx, method = "euclidean") # distance matrix
    fit <- hclust(d, method = "ward.D")
    cl.label<- cutree(fit, k = k) # cut tree into k clusters
  }

  if(method == 'mclust'){

    fit <- Mclust(reduced.mtx, G = k)
    cl.label = fit$classification
  }
  return(cl.label)
}



## do clustering using different #pcs and different methods, given an seurat.obj
## return seurat object, with a new meta.data column with name clusterLabelName
## gcSNN -- the clustering method used by seurat (Louvain algorithm on snn)
clust4seurat <- function(seurat.obj, npc = 30, method = 'gcSNN', reduction = 'tsne', resolution = 0.1, 
                                clustLabelName = paste0('clusterBy_', reduction, '_',  method), k = NULL){

  ## note the seurat object was done pca using 100 pcs; so npc should be smaller than 100
  ## using gcSNN
  if(is.null(seurat.obj@reductions$pca)){
    seurat.obj <- RunPCA(seurat.obj, npcs = npc, verbose = F)
  } 


  if(reduction == 'tsne'){
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc)
  }
  if(reduction == 'umap'){
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:npc, verbose = F)
  }



  if(method == 'gcSNN'){
    
    if(!is.null(k)){
      # find best resolution to get k cluster
      resolution = queryResolution4Seurat(seurat.obj, reduction = reduction, npc = npc, min_resl = resolution,
                                          max_resl = 10*resolution, k = k)
    }

    if(reduction == 'tsne'){
        
        seurat.obj <- FindNeighbors(seurat.obj, reduction = reduction, dims = 1:2)
        seurat.obj <- FindClusters(seurat.obj, resolution = resolution, verbose = F)
    }

    if(reduction == 'umap'){
        seurat.obj <- FindNeighbors(seurat.obj, reduction = reduction, dims = 1:2)
        seurat.obj <- FindClusters(seurat.obj, resolution = resolution, verbose = F)
    }

    if(reduction == 'pca'){
        seurat.obj <- FindNeighbors(seurat.obj, reduction = reduction, dims = 1:npc)
        seurat.obj <- FindClusters(seurat.obj, resolution = resolution, verbose = F)
    }
    seurat.obj@meta.data[[clustLabelName]] = as.integer(seurat.obj@active.ident)

  }else{
    if(reduction == 'tsne'){
       
        embd.mat <- seurat.obj@reductions$tsne@cell.embeddings
        
    }

    if(reduction == 'umap'){
        
        embd.mat <- seurat.obj@reductions$umap@cell.embeddings
    }

    if(reduction == 'pca'){
        embd.mat <- seurat.obj@reductions$pca@cell.embeddings[, 1:npc]
    }

     
     seurat.obj@meta.data[[clustLabelName]] = basicCluster(embd.mat, method, k)
  }


 
  return(seurat.obj)
}


# fit cistopic model and using the predicted cell * topic probability matrix
# for clustering
run_cisTopic <- function(mtx, nCores = 4, frac_in_cell = 0.025){
  # prepare the right format of rownames
  rnames = data.table('region' = rownames(mtx))
  tmp = tidyr::separate(rnames, col = 'region', into = c('chr', 'start', 'end'))
  rnames = paste0(tmp$chr, ':', tmp$start, '-', tmp$end)
  rownames(mtx) = rnames

  mtx0 = 1 * (mtx > 0)
  rr = Matrix::rowMeans(mtx0)
  mtx = mtx[rr >= frac_in_cell, ]
  cisTopicObject <- createcisTopicObject(mtx, project.name='scATAC')
  rm(mtx, mtx0)
  cisTopicObject <- runModels(cisTopicObject, topic = c(10, 20, 30, 40, 50, 80, 100), seed = 987, nCores = nCores, 
    burnin = 120, iterations = 150, addModels = T)
  #cisTopicObject <- selectModel(cisTopicObject, keepBinarymatrix = F, keepModels = F)
  #cellassign <- t(modelMatSelection(cisTopicObject, 'cell', 'Probability'))
  return(cisTopicObject)
}




run_chromVAR <- function(mtx, genomeName = 'BSgenome.Hsapiens.UCSC.hg38',
                         ncore = 3){
  
  register(MulticoreParam(ncore))
  if(!require(genomeName, character.only = T)) BiocManager::install(genomeName)
  
  peaks = data.table('x' = rownames(mtx))
  peaks = tidyr::separate(peaks, col = 'x', into = c('chr', 'start', 'end'), sep = '-')
  peaks = GenomicRanges::makeGRangesFromDataFrame(peaks)
  
  frag.counts = SummarizedExperiment(assay = list(counts = mtx),
                                     rowRanges = peaks)
  frag.counts <- addGCBias(frag.counts, genome = genomeName)
  #motifs <- getJasparMotifs()
  library(chromVARmotifs) ## cisbp motif
 
  if(grepl(genomeName, pattern = 'hg')){
    motifs = human_pwms_v2
  }else{
    motifs = mouse_pwms_v2
  }
  
  motif_ix <- matchMotifs(motifs, frag.counts,
                          genome = genomeName)
  dev <- computeDeviations(object = frag.counts, 
                           annotations = motif_ix)
  bg <- getBackgroundPeaks(object = frag.counts)
  
  dev <- computeDeviations(object = frag.counts, annotations = motif_ix,
                           background_peaks = bg)
  
  #motif.zscore = dev@assays$data$z
  return(dev)
}


## do graph cluster snn on cistop object
gcSNN4CistopicObj <- function(cisTopic.obj, ntopic = 20, clust_by = 'topic',
                            resolution = 0.1, k = NULL){
  #note cisTopic.obj should include model with ntopic
  sele.cisTopic <- selectModel(cisTopic.obj, select = ntopic, 
                               keepBinaryMatrix = F, keepModels = F)
  sele.cisTopic <- runtSNE(sele.cisTopic, target='cell')
  sele.cisTopic <- runUmap(sele.cisTopic, target='cell')
  
  if(clust_by == 'topic'){
    embd.mtx = t(modelMatSelection(sele.cisTopic, 'cell', 'Probability'))
  }
  if(clust_by == 'tsne'){
    embd.mtx = sele.cisTopic@dr$cell$tSNE
  }
  if(clust_by == 'umap'){
    embd.mtx = sele.cisTopic@dr$cell$Umap
  }
  
  snn = FindNeighbors(dist(embd.mtx))
  names(snn) = c('nn', 'snn')
  # creat a seurat object to faciliat finding clusters using seurat
  tmp.seurat <- CreateSeuratObject(sele.cisTopic@count.matrix)
  tmp.seurat@graphs = snn
  if(!is.null(k)) {

    resolution = queryResolution4Topic(tmp.seurat, k = k, min_resl = resolution, max_resl = 10 * resolution)
  }

  tmp.seurat <- FindClusters(tmp.seurat, graph.name = 'snn', resolution = resolution)
  #tmp.seurat@reductions$tsne@cell.embeddings = sele.cisTopic@dr$cell$tSNE
  #tmp.seurat@reductions$umap@cell.embeddings = sele.cisTopic@dr$cell$Umap
  
  colnames(sele.cisTopic@dr$cell$tSNE) = c('tSNE_1', 'tSNE_2')
  colnames(sele.cisTopic@dr$cell$Umap) = c('UMAP_1', 'UMAP_2')
  
  return(list('tsne_coord' = sele.cisTopic@dr$cell$tSNE, 
              'umap_coord' = sele.cisTopic@dr$cell$Umap,
              'cluster_label' = as.integer(tmp.seurat@active.ident)))
}


##calculate rand index, given true labels and npc, and a seurat object
## return data.table recoords rand index for different conditions
calRand_gTrueLabelsAndNPC <- function(seurat.obj, true.labels, npc = 20, 
                                      clust_methods = c('gcSNN', 'hclust', 'kmeans'), 
                                      reductions = c('pca', 'tsne', 'umap'), resolution = 0.1){
   
   if(is.null(seurat.obj@reductions$pca)){
    seurat.obj <- RunPCA(seurat.obj, npcs = npc, verbose = F)
  } 

   seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc)
   seurat.obj <- RunUMAP(seurat.obj, dims = 1:npc, verbose = F)

   rands = matrix(0, length(clust_methods), length(reductions))

   for(i in 1:nrow(rands)){
    for (j in 1:ncol(rands)){
      seurat.obj <- clust4seurat(seurat.obj, npc = npc, method = clust_methods[i], 
                                  reduction = reductions[j], resolution = resolution, k = length(unique(true.labels)),
                                  clustLabelName = paste(clust_methods[i], '_', reductions[j], '_npc', npc))
      rands[i, j] = adjustedRandIndex(seurat.obj@meta.data[[paste(clust_methods[i], '_', reductions[j], '_npc', npc)]], true.labels)
    }
   }
   rownames(rands) = clust_methods
   colnames(rands) = reductions
   rands = reshape2::melt(rands, value.name = 'rand')
   rands$npc = npc
   return(rands)
} 


##calculate rand index, given true labels and npc, and a seurat object
## return data.table recoords rand index for different conditions
calRand_gTrueLabelsAndNTOPIC <- function(cisTopic.obj, true.labels, ntopic = 20, 
                                      clust_methods = c('gcSNN', 'hclust', 'kmeans'), 
                                      reductions = c('topic', 'tsne', 'umap'), resolution = 0.1){
   

   rands = matrix(0, length(clust_methods), length(reductions))

   for(i in 1:nrow(rands)){

    if(clust_methods[i] == 'gcSNN'){
      for (j in 1:ncol(rands)){
      pred.labels <- gcSNN4CistopicObj(cisTopic.obj, ntopic = ntopic, clust_by = reductions[j], resolution = resolution, 
        k = length(unique(true.labels)))$cluster_label

      rands[i, j] = adjustedRandIndex(pred.labels, true.labels)
      }
    }else{
      for (j in 1:ncol(rands)){
        sele.cisTopic <- selectModel(cisTopic.obj, select = ntopic, 
                               keepBinaryMatrix = F, keepModels = F)
        sele.cisTopic <- runtSNE(sele.cisTopic, target='cell')
        sele.cisTopic <- runUmap(sele.cisTopic, target='cell')
        if(reductions[j] == 'topic') cell_topic <- t(modelMatSelection(sele.cisTopic, 'cell', 'Probability'))
        if(reductions[j] == 'tsne') cell_topic <- sele.cisTopic@dr$cell$tSNE
        if(reductions[j] == 'umap') cell_topic <- sele.cisTopic@dr$cell$Umap
        pred.labels <- basicCluster(cell_topic, method = clust_methods[i], k = length(unique(true.labels)))

        rands[i, j] = adjustedRandIndex(pred.labels, true.labels)
     }
   }
 }
   rownames(rands) = clust_methods
   colnames(rands) = reductions
   rands = reshape2::melt(rands, value.name = 'rand')
   rands$ntopic = ntopic

   return(rands)
} 


## compare clustering rand index using different npc, given a method and reduction
## plot tsne/umap or no plot (plotDR)
compRand_npc_gReductionAndMethod <- function(seurat.obj, clust_method = 'gcSNN', reduction = 'pca',
                                resolution = 0.2, k = 10, npcs = c(20, 30, 50, 100), plotDR = TRUE, ...){

  
  if(is.null(seurat.obj@reductions$pca)){
    seurat.obj <- RunPCA(seurat.obj, npcs = max(npcs), verbose = F)
  } else if(ncol(seurat.obj@reductions$pca@cell.embeddings) < max(npcs)){
    seurat.obj <- RunPCA(seurat.obj, npcs = max(npcs), verbose = F)
  }
  

  i = 0
  pp_tsne = pp_umap = list()
  for(npc0 in npcs){
    i = i + 1
    seurat.obj <- clust4seurat(seurat.obj, npc = npc0, method = clust_method, 
                                  reduction = reduction, resolution = resolution, k = k,
                                  clustLabelName = paste(clust_method, '_', reduction, '_npc', npc0))
    if(plotDR) {
       seurat.obj = RunTSNE(seurat.obj, dims = 1:npc0)
      pp_tsne[[i]] = DimPlot(seurat.obj, reduction = 'tsne', group.by = paste(clust_method, '_', reduction, '_npc', npc0))
       seurat.obj = RunUMAP(seurat.obj, dims = 1:npc0, verbose = F)
      pp_umap[[i]] = DimPlot(seurat.obj, reduction = 'umap', group.by = paste(clust_method, '_', reduction, '_npc', npc0))

    }

  }
  #embd.mtx = seurat.obj@reductions[[reduction]]@cell.embeddings

  # evaluat clustering, calculate corrected rand index
  rands = rep(0, length(npcs) - 1)
  for(i in 1:length(rands)){
    #rands[i] = eval_cluster4mat(embd.mtx, 
    #seurat.obj@meta.data[[paste(clust_method, '_', reduction, '_npc', npcs[i])]],
    #seurat.obj@meta.data[[paste(clust_method, '_', reduction, '_npc', npcs[i +1])]])$rand
    rands[i] = adjustedRandIndex(seurat.obj@meta.data[[paste(clust_method, '_', reduction, '_npc', npcs[i])]],
    seurat.obj@meta.data[[paste(clust_method, '_', reduction, '_npc', npcs[i +1])]])
  }
  
  set.cols = brewer.pal(n = length(rands), name = 'Dark2')
  
  
  res = list('rand' = rands)
  if(plotDR){
    res$plots_tsne = pp_tsne
    res$plots_umap = pp_umap
  }

  return(res)

}


## compare clustering using different reductions, given npc and method  -- not used
compRand_reduction_gNPCAndMethod <- function(seurat.obj, clust_method = 'gcSNN', 
                                  reductions = c('pca', 'tsne', 'umap'),
                                  resolution = 0.2, k = 10, npc = 20, ...){

  if(is.null(seurat.obj@reductions$pca)){
    seurat.obj <- RunPCA(seurat.obj, npcs = npc, verbose = F)
  }

  for(reduction0 in reductions){
    seurat.obj <- clust4seurat(seurat.obj, npc = npc, method = clust_method, 
                                  reduction = reduction0, resolution = resolution, k = k,
                                  clustLabelName = paste(clust_method, '_', reduction0, '_npc', npc))

  }
  #embd.mtx = seurat.obj@reductions[['tsne']]@cell.embeddings

  # evaluat clustering, calculate corrected rand index
  rands = rep(0, length(reductions) - 1)
  for(i in 1:(length(rands))){
    #rands[i] = eval_cluster4mat(embd.mtx, 
    #seurat.obj@meta.data[[paste(clust_method, '_', reductions[i], '_npc', npc)]],
    #seurat.obj@meta.data[[paste(clust_method, '_', reductions[i+1], '_npc', npc)]])$rand
    rands[i] = adjustedRandIndex(seurat.obj@meta.data[[paste(clust_method, '_', reductions[i], '_npc', npc)]],
    seurat.obj@meta.data[[paste(clust_method, '_', reductions[i+1], '_npc', npc)]])
  }
  
  set.cols = brewer.pal(n = length(rands), name = 'Dark2')
  
  p <- barplot(rands, col = set.cols, ylab = 'Adjust Rand Index', ylim = c(0, 1),  ...)
  
  return(rands)
   
}

## compare clustering using different methods, given npc and reduction 
compRand_method_gNPCAndReduction <- function(seurat.obj, clust_methods = c('gcSNN', 'hclust', 'kmean'), 
                                  reduction = 'pca',
                                  resolution = 0.2, k = 10, npc = 20, plotDR = FALSE, ...){
  if(is.null(seurat.obj@reductions$pca)){
    seurat.obj <- RunPCA(seurat.obj, npcs = npc, verbose = F)
  } 
  i = 0
  pp_tsne = pp_umap = list()
  for(method0 in clust_methods){
    i = i + 1
    seurat.obj <- clust4seurat(seurat.obj, npc = npc, method = method0, 
                                  reduction = reduction, resolution = resolution, k = k,
                                  clustLabelName = paste(method0, '_', reduction, '_npc', npc))
    if(plotDR) {
       seurat.obj = RunTSNE(seurat.obj, dims = 1:npc)
      pp_tsne[[i]] = DimPlot(seurat.obj, reduction = 'tsne', group.by = paste(method0, '_', reduction, '_npc', npc))
       seurat.obj = RunUMAP(seurat.obj, dims = 1:npc, verbose = F)
      pp_umap[[i]] = DimPlot(seurat.obj, reduction = 'umap', group.by = paste(method0, '_', reduction, '_npc', npc))

    }
  }
  #embd.mtx = seurat.obj@reductions[[reduction]]@cell.embeddings

  # evaluat clustering, calculate corrected rand index
  rands = rep(0, length(clust_methods) - 1)
  for(i in 1:(length(rands))){
    #rands[i] = eval_cluster4mat(embd.mtx, 
    #seurat.obj@meta.data[[paste(clust_methods[i], '_', reduction, '_npc', npc)]],
    #seurat.obj@meta.data[[paste(clust_methods[i + 1], '_', reduction, '_npc', npc)]])$rand
    rands[i] = adjustedRandIndex(seurat.obj@meta.data[[paste(clust_methods[i], '_', reduction, '_npc', npc)]],
    seurat.obj@meta.data[[paste(clust_methods[i + 1], '_', reduction, '_npc', npc)]])

  }
  
  res = list('rand' = rands)
  if(plotDR){
    res$plots_tsne = pp_tsne
    res$plots_umap = pp_umap
  }
  return(res)

}




## compare clustering rand index using different npc, given a method and reduction
compRand_ntopic_gReductionAndMethod <- function(cisTopic.obj, clust_method = 'gcSNN', reduction = 'topic',
                                resolution = 0.2, k = 10, ntopics = c(20, 30, 40, 50), ...){

  
  snn.res = list()
  i=1
  for(ntopic0 in ntopics){
    if(clust_method == 'gcSNN') {
      snn.res[[i]] <- gcSNN4CistopicObj(cisTopic.obj, ntopic0, reduction, resolution = resolution, k = k)
      }else{
        
        

        sele.cisTopic <- selectModel(cisTopic.obj, select = ntopic0, 
                               keepBinaryMatrix = F, keepModels = F)
        sele.cisTopic <- runtSNE(sele.cisTopic, target='cell')
        sele.cisTopic <- runUmap(sele.cisTopic, target='cell')
        if(reduction == 'topic') cell_topic <- t(modelMatSelection(sele.cisTopic, 'cell', 'Probability'))
        if(reduction == 'tsne') cell_topic <- sele.cisTopic@dr$cell$tSNE
        if(reduction == 'umap') cell_topic <- sele.cisTopic@dr$cell$Umap
        snn.res[[i]] <- basicCluster(cell_topic, method = clust_method, k = k)

      }
    
    i = i+1
    
  }

  # evaluat clustering, calculate corrected rand index
  rands = rep(0, length(ntopics) - 1)
  for(i in 1:length(rands)){
    if(clust_method == 'gcSNN') {
        rands[i] = adjustedRandIndex(snn.res[[i]]$cluster_label,
                                snn.res[[i+1]]$cluster_label)
    }else{
        rands[i] = adjustedRandIndex(snn.res[[i]],
                                snn.res[[i+1]])
    }
  }
  
  set.cols = brewer.pal(n = length(rands), name = 'Dark2')
  
  p <- barplot(rands, col = set.cols, ylab = 'Adjust Rand Index', ylim = c(0, 1), ...)
  
  return(rands)

}


## compare clustering rand index using different npc, given a method and reduction
compRand_method_gNTOPICAndReduction <- function(cisTopic.obj, clust_methods = c('gcSNN', 'hclust', 'kmeans'), reduction = 'topic',
                                resolution = 0.2, k = 10, ntopic = 30, ...){

  
  snn.res = list()
  i=1
  for(method0 in clust_methods){
    if(clust_methods[i] == 'gcSNN') {
      snn.res[[i]] <- gcSNN4CistopicObj(cisTopic.obj, ntopic, reduction, resolution = resolution, k = k)
      }else{
        
        

        sele.cisTopic <- selectModel(cisTopic.obj, select = ntopic, 
                               keepBinaryMatrix = F, keepModels = F)
        sele.cisTopic <- runtSNE(sele.cisTopic, target='cell')
        sele.cisTopic <- runUmap(sele.cisTopic, target='cell')
        if(reduction == 'topic') cell_topic <- t(modelMatSelection(sele.cisTopic, 'cell', 'Probability'))
        if(reduction == 'tsne') cell_topic <- sele.cisTopic@dr$cell$tSNE
        if(reduction == 'umap') cell_topic <- sele.cisTopic@dr$cell$Umap
        snn.res[[i]] <- basicCluster(cell_topic, method = method0, k = k)

      }
    
    i = i+1
    
  }


  # evaluat clustering, calculate corrected rand index
  rands = rep(0, length(clust_methods) - 1)
  for(i in 1:length(rands)){
    if(clust_methods[i] == 'gcSNN') {
        rands[i] = adjustedRandIndex(snn.res[[i]]$cluster_label,
                                snn.res[[i+1]]$cluster_label)
    }else{
        rands[i] = adjustedRandIndex(snn.res[[i]],
                                snn.res[[i+1]])
    }
  }
  
  set.cols = brewer.pal(n = length(rands), name = 'Dark2')
  
  p <- barplot(rands, col = set.cols, ylab = 'Adjust Rand Index', ylim = c(0, 1), ...)
  
  return(rands)

}


## plot cisTopic cluster on dimension reduction plot
## transfer to seurat object to have a similar color scheme
plotDR4cisTopic_bySeurat <- function(cisTopic.obj, ntopic = 20, cistopic.cl.label){
  sele.cisTopic <- selectModel(cisTopic.obj, select = ntopic, 
                               keepBinaryMatrix = F, keepModels = F)
  sele.cisTopic <- runtSNE(sele.cisTopic, target='cell')
  sele.cisTopic <- runUmap(sele.cisTopic, target='cell')
  colnames(sele.cisTopic@dr$cell$tSNE) = c('tSNE_1', 'tSNE_2')
  colnames(sele.cisTopic@dr$cell$Umap) = c('UMAP_1', 'UMAP_2')
  


  embd.mtx = modelMatSelection(sele.cisTopic, 'cell', 'Probability')

  seurat.obj <- CreateSeuratObject(embd.mtx)
  seurat.obj@meta.data[['cistopic.cl.label']] = cistopic.cl.label
  seurat.obj = ScaleData(seurat.obj, do.scale = F, do.center = F, vars.to.regress = NULL)
  seurat.obj = RunPCA(seurat.obj, npc = floor(ntopic/2), features = rownames(seurat.obj), verbose = F)


  seurat.obj = RunTSNE(seurat.obj, dims = 1:floor(ntopic/2))
  #replace the tsne coordinate
  seurat.obj@reductions$tsne@cell.embeddings = sele.cisTopic@dr$cell$tSNE
  pp_tsne = DimPlot(seurat.obj, reduction = 'tsne', group.by = 'cistopic.cl.label')

  
  seurat.obj = RunUMAP(seurat.obj, dims = 1:floor(ntopic/2), verbose = F)
  seurat.obj@reductions$umap@cell.embeddings = sele.cisTopic@dr$cell$Umap
  #replace the umap coordinate
  pp_umap = DimPlot(seurat.obj, reduction = 'umap', group.by = 'cistopic.cl.label')



  return(list('plot_tsne' = pp_tsne, 'plot_umap' = pp_umap, 'seurat.obj' = seurat.obj))
}

## get overlap number of nearby DAs of DE within a given dist
getOverlap_DEwithDA <- function(DEs, DAs, flank.dist = 500 * 1000){
  chrs = unique(DEs$chr)
  res.overlap = NULL
  for(chr0 in chrs){
    DE0 = DEs[chr == chr0, ]
    DA0 = DAs[chr == chr0, ]
    if(nrow(DA0) == 0) next
    tmp.res = lapply(DE0$midP, function(x) DA0[abs(midP - x) <= flank.dist])
    tmp.res = do.call('rbind', tmp.res)
    res.overlap = rbind(res.overlap, tmp.res)
  }
  res.overlap = res.overlap[!duplicated(res.overlap), ]
  return(res.overlap)
}

## link DEs with DAs 
## given a table of DE and a tables DAs (with cluster info), output a matrix records the number of 
## DAs for each groups DEs with a flank.dist 
## the DE_table should include columns <cluster> <gene> <p_val>
## and DA_table should include columns <cluster> <peak> <p_val>
## Note: peak should be like chr-start-end, p_val or logFC should be 
linkDEwithDA <- function(DE_table, DA_table, gene_ann, flank.dist = 500*1000, 
                      top.de = 0.2, top.da = 0.2, sort_by = 'p_val', enhs = NULL){
  de.cl = unique(DE_table$cluster)
  da.cl = unique(DA_table$cluster)
  len1 = length(de.cl)
  len2 = length(da.cl)
  noverlaps = noverlaps.a2e = matrix(0, len1, len2)
 
  # add coordinate to DE
  gene_ann = gene_ann[!duplicated(gene_name), ]
  setkey(gene_ann, gene_name)
  DE_table = DE_table[gene %in% gene_ann$gene_name, ]
  DE_table[, 'chr' := gene_ann[J(DE_table$gene)]$chr]
  DE_table[, 'midP' := gene_ann[J(DE_table$gene)]$tss]
  
  # get midpoint of DA
  DA_table[, 'chr':= unlist(strsplit(peak, '-'))[1], by = peak]
  DA_table[, 'start':= as.integer(unlist(strsplit(peak, '-'))[2]), by = peak]
  DA_table[, 'end':= as.integer(unlist(strsplit(peak, '-'))[3]), by = peak]
  DA_table[, 'midP' := floor(start/2 + end/2)]
  
  # only use peak overlapped with enhancers
  if(!is.null(enhs)){
    filter.da = NULL
    chrs = unique(DA_table$chr)
    for(chr0 in chrs){
      enh0 = enhs[chr == chr0, ]
      DA0 = DA_table[chr == chr0, ]
      DA0[, 'mdist' := min(abs(midP - enh0$midP)), by = midP]
      DA0[, 'len' := end - start]
      
      filter.da = rbind(filter.da, DA0[mdist < 1000 + len/2])
    }
    filter.da[, c('mdist', 'len') := NULL]
    DA_table = filter.da
    rm(filter.da)
  }
  
  # calculate overlaps
  for(i in 1:len1){
    DEs = DE_table[cluster == de.cl[i], ]
    setkeyv(DEs, sort_by)
    for(j in 1:len2){
      DAs = DA_table[cluster == da.cl[j], ]
      setkeyv(DAs, sort_by)
      topDE = floor(nrow(DEs) * top.de)
      topDA = floor(nrow(DAs) * top.da)
      noverlaps[i, j] = nrow(getOverlap_DEwithDA(DEs[1:topDE, ], DAs[1:topDA, ],
                                            flank.dist))
      noverlaps.a2e[i, j] = nrow(getOverlap_DEwithDA(DAs[1:topDA, ], DEs[1:topDE, ], 
                                                flank.dist))
    }
  }
  
  # transform to pvalue(hyper-geometric test), given overlap mat and
  # col and row sums (union)
  get_pv <- function(noverlaps, da.count, de.count, nn){
    #da.count = colSums(noverlaps)
    #de.count = rowSums(noverlaps)
    
    da.mat = matrix(rep(da.count, each = nrow(noverlaps)), 
                    nrow = nrow(noverlaps))
    de.mat = matrix(rep(de.count, each = ncol(noverlaps)), 
                    nrow = nrow(noverlaps), byrow = T)
    
    pvs = phyper(noverlaps, da.mat, nn - da.mat,
                 de.mat, lower.tail = F)
    # normalized overlaps
    noverlaps.normDA = noverlaps/(da.mat * de.mat) * nn
    #tt = t(percentize(t(noverlaps.normDA)))
    
    return(list(pvs, noverlaps.normDA))
  }
  
 
  col.pool = lapply(da.cl, function(x) getOverlap_DEwithDA(DE_table, DA_table[cluster == x, ],  flank.dist))
  row.pool = lapply(de.cl, function(x) getOverlap_DEwithDA(DE_table[cluster == x, ], DA_table,  flank.dist))
  all.over = do.call('rbind', col.pool)
  all.over = all.over[!duplicated(all.over), ] ## should remove duplicates
  nn = nrow(all.over) # total 'balls'
  
  res = get_pv(noverlaps, sapply(col.pool, function(x) nrow(x)), 
               sapply(row.pool, function(x) nrow(x)), nn)
  
  col.pool = lapply(da.cl, function(x) getOverlap_DEwithDA(DA_table[cluster == x, ], DE_table, flank.dist))
  row.pool = lapply(de.cl, function(x) getOverlap_DEwithDA(DA_table, DE_table[cluster == x, ],  flank.dist))
  all.over = do.call('rbind', col.pool)
  all.over = all.over[!duplicated(all.over), ]
  nn = nrow(all.over)
  res.a2e = get_pv(noverlaps.a2e, sapply(col.pool, function(x) nrow(x)), 
                   sapply(row.pool, function(x) nrow(x)), nn)
  
  return(list('pvs.e2a' = res[[1]], 'norm.e2a' = res[[2]], 
               'pvs.a2e' = res.a2e[[1]], 'norm.a2e' = res.a2e[[2]],
              'count.e2a' = noverlaps, 'count.a2e' = noverlaps.a2e,
              'DA_filtered' = DA_table))
  
}


# mtx: gene by cell matrix, or peak/bin by cell matrix
# seurat.obj: optional, if not provided, will construct one (which will clustering on pca1:50),
# if seurat.obj is provided, assume it did pca
# norm_mtx: normalized matrix, equals to log(mtx +1) if not provided
# extraDims: do extra dimension reduction on 10, 30, 100  etc
# subSamples: provided a subsample version
# vFeatures: given variable features; if NULL using default variable features
prepInput4Cello <- function(mtx, seurat.obj, norm_mtx = NULL, 
                            cello.name = 'scRNA', assay = 'RNA', 
                            resl4clust = 0.6, 
                            extraDims = c(10, 20, 30, 50, 80, 100),
                            subSamples = NULL, subCelloName = 'sub',
                            vars.to.regOnPca = NULL, downSample = NULL,
                            vFeatures = NULL){

    if(!is.null(vFeatures)) seurat.obj <- ScaleData(seurat.obj, features = vFeatures)
    if(is.null(vFeatures)) vFeatures = VariableFeatures(seurat.obj)  
    DefaultAssay(seurat.obj) = assay
    
    ndefault = ncol(seurat.obj@reductions$pca@cell.embeddings)
    
    seurat.obj <- RunPCA(seurat.obj, npcs = ndefault, verbose = F,
                         assay = assay, features = vFeatures)
    
    if(!is.null(vars.to.regOnPca)) seurat.obj = regress_on_pca(seurat.obj, vars.to.regOnPca)
    
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:ndefault, check_duplicates = FALSE,
                          assay = assay)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:ndefault, verbose = F,
                          assay = assay)
    my_tsne_proj <- seurat.obj@reductions$tsne@cell.embeddings
    my_umap_proj <- seurat.obj@reductions$umap@cell.embeddings 
    
    
    
    if(ndefault < 100){
      
      seurat.obj <- RunPCA(seurat.obj, npcs = 100, verbose = F,
                           assay = assay, features = vFeatures)
      
      if(!is.null(vars.to.regOnPca)) seurat.obj = regress_on_pca(seurat.obj, vars.to.regOnPca)
      
      seurat.obj <- RunTSNE(seurat.obj, dims = 1:100, check_duplicates = FALSE,
                            assay = assay)
      seurat.obj <- RunUMAP(seurat.obj, dims = 1:100, verbose = F,
                            assay = assay)
      my_tsne_proj100 <- seurat.obj@reductions$tsne@cell.embeddings
      my_umap_proj100 <- seurat.obj@reductions$umap@cell.embeddings 
      
      
      
    }
  
  
  
  #mtx = mtx[rownames(mtx) %in% rownames(seurat.obj), ]
  mtx = mtx[, colnames(mtx) %in% colnames(seurat.obj)]
  
  meta = seurat.obj@meta.data
  fmeta <- data.frame(symbol = rownames(mtx)) 
  rownames(fmeta) <- fmeta$symbol
  
  # preparing expression matrix
  if(is.null(norm_mtx)) norm_mtx = log2(mtx + 1)
  
  ids = 1:ncol(seurat.obj)
  if(!is.null(downSample)){
    set.seed(2019)
    if(downSample <= ncol(seurat.obj)) 
      ids = sample(1:ncol(seurat.obj), downSample)
  }
  
  eset <- new("ExpressionSet",
              assayData = assayDataNew("environment", 
                                       exprs = mtx[, ids], 
                                       norm_exprs = norm_mtx[, ids]),
              phenoData =  new("AnnotatedDataFrame", data = meta[ids, ]),
              featureData = new("AnnotatedDataFrame", data = fmeta))
  
  # Creating a cello 
  Cello <- setClass("Cello",
                    slots = c(
                      name = "character", # The name of cvis
                      idx = "numeric", # The index of global cds object
                      proj = "list", # The projections as a list of data frame
                      pmeta = "data.frame", # The local meta data
                      notes = "character" # Other information to display to the user
                    )
  )
  
  cello1 <- new("Cello", name = cello.name, idx = 1:length(ids)) 
  
  my_pca_proj100 <- seurat.obj@reductions$pca@cell.embeddings
  my_tsne_proj100 <- seurat.obj@reductions$tsne@cell.embeddings
  my_umap_proj100 <- seurat.obj@reductions$umap@cell.embeddings
  
  cello1@proj <- list("PCA" = my_pca_proj100[ids, ]) 
  if(is.null(extraDims)) extraDims = unique(c(ndefault, 100))
  
  for(dim0 in extraDims){
    if(dim0 != 100 & dim0 != ndefault){
    
      seurat.obj <- RunTSNE(seurat.obj, dims = 1:dim0, check_duplicates = FALSE,
                          assay = assay)
      seurat.obj <- RunUMAP(seurat.obj, dims = 1:dim0, verbose = F,
                          assay = assay)
      my_tsne_proj0 <- seurat.obj@reductions$tsne@cell.embeddings
      my_umap_proj0 <- seurat.obj@reductions$umap@cell.embeddings
    }
    
    if(dim0 == ndefault){
      my_tsne_proj0 <- my_tsne_proj
      my_umap_proj0 <- my_umap_proj
    }
    
    if(dim0 == 100 & dim0 != ndefault){
      my_tsne_proj0 <- my_tsne_proj100
      my_umap_proj0 <- my_umap_proj100
    }
    
    cello1@proj[[paste0('t-SEN', dim0)]] <- my_tsne_proj0[ids, ]
    cello1@proj[[paste0('UMAP', dim0)]] <- my_umap_proj0[ids, ]
   
  }
 
  #cello1@proj = cello1@proj[order(names(cello1@proj))]

  clist <- list()
  clist[[cello.name]] <- cello1
 
  if(!is.null(subSamples)){
    ids = which(colnames(seurat.obj) %in% subSamples)
    cello2 <- new("Cello", name = cello.name, idx = ids) 
    
    
    cello2@proj <- list("PCA" = my_pca_proj[ids, ], 
                        't-SNE50' = my_tsne_proj50[ids, ],
                        "UMAP50" = my_umap_proj50[ids, ]) 
   
    cello.name.sub = paste0(cello.name, '_', subCelloName)
    clist[[cello.name.sub]] <- cello2
  }
  
  
  return(list('eset' = eset, 'clist' = clist))
}

## update different version of clustering, and other information
## suppose cells are the same in seurat.obj, and only 1 clist obj
## updateProjDims means update tsne and umap using pca updateProjDims
## add extra meta data for each cell
updateInput4Cello <- function(celloInputPath, seuratObjPath, 
                              assay = 'integrated', extraMeta = NULL,
                              addPcaDims = c(20),
                              clusterOn = 'pca', clusterOnDim = 30,
                              resolutions = c(0.2, 0.4, 0.6)){
  eset = readRDS(paste0(celloInputPath, '/eset.rds'))
  clist = readRDS(paste0(celloInputPath, '/clist.rds'))
  mdata = phenoData(eset)
  
  seurat.obj = readRDS(seuratObjPath)
  DefaultAssay(seurat.obj) <- assay
  
  
  
  ## check whether clusterOnDim exists
  if(!is.null(clusterOnDim)){
    if(ncol(seurat.obj@reductions$pca@cell.embeddings) < clusterOnDim){
      seurat.obj <- RunPCA(seurat.obj, npcs = clusterOnDim, verbose = F)
      seurat.obj <- RunTSNE(seurat.obj, dims = 1:clusterOnDim, check_duplicates = FALSE,
                            assay = assay)
      seurat.obj <- RunUMAP(seurat.obj, dims = 1:clusterOnDim, verbose = F,
                            assay = assay)
      my_tsne_proj <- seurat.obj@reductions$tsne@cell.embeddings
      my_umap_proj <- seurat.obj@reductions$umap@cell.embeddings 
      
      clist[[1]]@proj[[paste0('t-SNE', clusterOnDim)]] = my_tsne_proj
      clist[[1]]@proj[[paste0('UMAP', clusterOnDim)]] = my_umap_proj
      
    }
    
  }
  
  ## try different clustering
  if(!is.null(clusterOnDim)){
    seurat.obj <- FindNeighbors(seurat.obj, reduction = clusterOn, dims = 1:clusterOnDim)
    for(res0 in resolutions){
      seurat.obj <- FindClusters(seurat.obj, resolution = res0)
      mdata@data[[paste0(assay, '_snn_', clusterOn, clusterOnDim , 'res.', res0)]] <- 
        seurat.obj[[paste0(assay, '_snn_res.', res0)]][, 1]
    }
  }
  
  mdata@data$seurat_clusters = seurat.obj$seurat_clusters
  
  if(!is.null(extraMeta)){
    if(all(rownames(extraMeta) == rownames(mdata))) {
      shared.features = intersect(colnames(mdata@data), colnames(extraMeta))
      mdata@data[, shared.features] <- NULL
      mdata@data = cbind(mdata@data, extraMeta)
    }
  }
  
  phenoData(eset) <- mdata
  
  
  ## update clist
  for(dim0 in addPcaDims){
    if(ncol(seurat.obj@reductions$pca@cell.embeddings) < dim0) {
      seurat.obj <- RunPCA(seurat.obj, verbose = F, npcs = dim0)
    }
    
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:dim0, check_duplicates = FALSE,
                          assay = assay)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:dim0, verbose = F,
                          assay = assay)
    my_tsne_proj <- seurat.obj@reductions$tsne@cell.embeddings
    my_umap_proj <- seurat.obj@reductions$umap@cell.embeddings 
    
    clist[[1]]@proj[[paste0('t-SNE', dim0)]] = my_tsne_proj
    clist[[1]]@proj[[paste0('UMAP', dim0)]] = my_umap_proj
  }
  
  
  saveRDS(eset, paste0(celloInputPath, '/eset.rds'))
  saveRDS(clist, paste0(celloInputPath, '/clist.rds'))
  if(!is.null(clusterOnDim) & !is.null(addPcaDims)) saveRDS(seurat.obj, seuratObjPath)
}



## update mtx of cello objects, without touch other stuff
## need provide feature names
updateRow4Cello <- function(celloInputPath, seuratObjPath, features_name,
                              assay = 'ATAC'){
  eset = readRDS(paste0(celloInputPath, '/eset.rds'))
  mdata = phenoData(eset)
  
  seurat.obj = readRDS(seuratObjPath)
  DefaultAssay(seurat.obj) <- assay
  
  mtx = seurat.obj@assys[[assay]]@counts[features_name, ]
  norm_mtx = seurat.obj@assys[[assay]]@data[features_name, ]
  
  meta = seurat.obj@meta.data
  fmeta <- data.frame(symbol = rownames(mtx)) 
  rownames(fmeta) <- fmeta$symbol
  
  # preparing expression matrix
  eset <- new("ExpressionSet",
              assayData = assayDataNew("environment", 
                                       exprs = mtx, 
                                       norm_exprs = norm_mtx),
              phenoData =  new("AnnotatedDataFrame", data = meta),
              featureData = new("AnnotatedDataFrame", data = fmeta))
  
  
  phenoData(eset) <- mdata
  
  saveRDS(eset, paste0(celloInputPath, '/eset.rds'))
 }



prepInput4Cello_harmony <- function(mtx, seurat.obj, norm_mtx = NULL, 
                            cello.name = 'scRNA', assay = 'RNA', 
                            resl4clust = 0.6, 
                            extraDims = c(10, 20, 30, 50, 80, 100),
                            subSamples = NULL, subCelloName = 'sub',
                            vars.to.harmony = NULL, downSample = NULL,
                            vFeatures = NULL){
  
  if(!is.null(vFeatures)) seurat.obj <- ScaleData(seurat.obj, features = vFeatures)
  if(is.null(vFeatures)) vFeatures = VariableFeatures(seurat.obj)  
  DefaultAssay(seurat.obj) = assay
  
  ndefault = ncol(seurat.obj@reductions$pca@cell.embeddings)
  
  seurat.obj <- RunPCA(seurat.obj, npcs = ndefault, verbose = F,
                       assay = assay, features = vFeatures)
  
  seurat.obj <- RunHarmony(seurat.obj, vars.to.harmony)
  seurat.obj <- RunTSNE(seurat.obj, dims = 1:ndefault, check_duplicates = FALSE,
                        assay = assay, reduction = 'harmony')
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:ndefault, verbose = F,
                        assay = assay, reduction = 'harmony')
  my_tsne_proj <- seurat.obj@reductions$tsne@cell.embeddings
  my_umap_proj <- seurat.obj@reductions$umap@cell.embeddings 
  
  
  
  if(ndefault < 100){
    
    seurat.obj <- RunPCA(seurat.obj, npcs = 100, verbose = F,
                         assay = assay, features = vFeatures)
    seurat.obj <- RunHarmony(seurat.obj, vars.to.harmony)
     
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:100, check_duplicates = FALSE,
                          assay = assay, reduction = 'harmony')
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:100, verbose = F,
                          assay = assay, reduction = 'harmony')
    my_tsne_proj100 <- seurat.obj@reductions$tsne@cell.embeddings
    my_umap_proj100 <- seurat.obj@reductions$umap@cell.embeddings 
    
    
    
  }
  
  
  
  #mtx = mtx[rownames(mtx) %in% rownames(seurat.obj), ]
  mtx = mtx[, colnames(mtx) %in% colnames(seurat.obj)]
  
  meta = seurat.obj@meta.data
  fmeta <- data.frame(symbol = rownames(mtx)) 
  rownames(fmeta) <- fmeta$symbol
  
  # preparing expression matrix
  if(is.null(norm_mtx)) norm_mtx = log2(mtx + 1)
  
  ids = 1:ncol(seurat.obj)
  if(!is.null(downSample)){
    set.seed(2019)
    if(downSample <= ncol(seurat.obj)) 
      ids = sample(1:ncol(seurat.obj), downSample)
  }
  
  eset <- new("ExpressionSet",
              assayData = assayDataNew("environment", 
                                       exprs = mtx[, ids], 
                                       norm_exprs = norm_mtx[, ids]),
              phenoData =  new("AnnotatedDataFrame", data = meta[ids, ]),
              featureData = new("AnnotatedDataFrame", data = fmeta))
  
  # Creating a cello 
  Cello <- setClass("Cello",
                    slots = c(
                      name = "character", # The name of cvis
                      idx = "numeric", # The index of global cds object
                      proj = "list", # The projections as a list of data frame
                      pmeta = "data.frame", # The local meta data
                      notes = "character" # Other information to display to the user
                    )
  )
  
  cello1 <- new("Cello", name = cello.name, idx = 1:length(ids)) 
  
  my_harmony_proj100 <- seurat.obj@reductions$harmony@cell.embeddings[, 1:2]
  my_tsne_proj100 <- seurat.obj@reductions$tsne@cell.embeddings
  my_umap_proj100 <- seurat.obj@reductions$umap@cell.embeddings
  
  cello1@proj <- list("PCA" = my_harmony_proj100[ids, ]) 
  if(is.null(extraDims)) extraDims = unique(c(ndefault, 100))
  
  for(dim0 in extraDims){
    if(dim0 != 100 & dim0 != ndefault){
      
      seurat.obj <- RunTSNE(seurat.obj, dims = 1:dim0, check_duplicates = FALSE,
                            assay = assay, reduction = 'harmony')
      seurat.obj <- RunUMAP(seurat.obj, dims = 1:dim0, verbose = F,
                            assay = assay, reduction = 'harmony')
      my_tsne_proj0 <- seurat.obj@reductions$tsne@cell.embeddings
      my_umap_proj0 <- seurat.obj@reductions$umap@cell.embeddings
    }
    
    if(dim0 == ndefault){
      my_tsne_proj0 <- my_tsne_proj
      my_umap_proj0 <- my_umap_proj
    }
    
    if(dim0 == 100 & dim0 != ndefault){
      my_tsne_proj0 <- my_tsne_proj100
      my_umap_proj0 <- my_umap_proj100
    }
    
    cello1@proj[[paste0('t-SEN', dim0)]] <- my_tsne_proj0[ids, ]
    cello1@proj[[paste0('UMAP', dim0)]] <- my_umap_proj0[ids, ]
    
  }
  
  #cello1@proj = cello1@proj[order(names(cello1@proj))]
  
  clist <- list()
  clist[[cello.name]] <- cello1
  
  if(!is.null(subSamples)){
    ids = which(colnames(seurat.obj) %in% subSamples)
    cello2 <- new("Cello", name = cello.name, idx = ids) 
    
    
    cello2@proj <- list("harmony" = my_harmony_proj[ids, ], 
                        't-SNE50' = my_tsne_proj50[ids, ],
                        "UMAP50" = my_umap_proj50[ids, ]) 
    
    cello.name.sub = paste0(cello.name, '_', subCelloName)
    clist[[cello.name.sub]] <- cello2
  }
  
  
  return(list('eset' = eset, 'clist' = clist))
}


# mtx: gene cis activity matrix (gene promoter + gene body)
# seurat.obj: should provided, if not provided, will construct one (which will clustering on pca1:50),
# if seurat.obj is provided
# norm_mtx: normalized matrix, equals to log(mtx +1) if not provided
# extraDR: do extra dimension reduction on 10, 30, 100 pcs; default 50
# cluster_excluded: provided a clean version with some clusters excluded
prepInput4Cello_gene_cis <- function(mtx, seurat.obj, 
                                          norm_mtx = NULL,
                            cello.name = 'scATAC_GA', assay = 'ATAC', 
                            resl4clust = 0.4,
                            top.variable = 0.1, extraDR = T,
                            cluster_excluded = NULL){
  
  if(is.null(seurat.obj)){
    #create a seurat object
    
    seurat.obj <- doBasicSeurat_new(mtx, doLog = T,  assay = assay,
                                  reg.var = paste0('nCount_ATAC'),
                                  top.variable = top.variable, npc = 30)
      
    
    
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:30, check_duplicates = FALSE,
                          assay = assay)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:30, verbose = F,
                          assay = assay)
    seurat.obj <- FindNeighbors(seurat.obj, reduction = 'pca', dims = 1:30,
                                assay = assay)
    seurat.obj <- FindClusters(seurat.obj, res = resl4clust)
    
  }else{
    DefaultAssay(seurat.obj) = assay
    if(is.null(seurat.obj@reductions$pca) || ncol(seurat.obj@reductions$pca@cell.embeddings) < 30){
      seurat.obj <- RunPCA(seurat.obj, npcs = 30, verbose = F,
                           assay = assay)
      seurat.obj <- RunTSNE(seurat.obj, dims = 1:30, check_duplicates = FALSE,
                            assay = assay)
      seurat.obj <- RunUMAP(seurat.obj, dims = 1:30, verbose = F,
                            assay = assay)
      if(is.null(seurat.obj@meta.data$seurat_clusters)){
        seurat.obj <- FindNeighbors(seurat.obj, reduction = 'pca', dims = 1:30,
                                    assay = assay)
        seurat.obj <- FindClusters(seurat.obj, res = resl4clust)
      }
      
    }
  }
  
  
  #mtx = mtx[rownames(mtx) %in% rownames(seurat.obj), ]
  mtx = mtx[, colnames(mtx) %in% colnames(seurat.obj)]
  
  meta = seurat.obj@meta.data
  fmeta <- data.frame(symbol = rownames(mtx)) 
  rownames(fmeta) <- fmeta$symbol
  
  # preparing expression matrix
  if(is.null(norm_mtx)) norm_mtx = log2(mtx + 1)
  eset <- new("ExpressionSet",
              assayData = assayDataNew("environment", 
                                       exprs = mtx, 
                                       norm_exprs = norm_mtx),
              phenoData =  new("AnnotatedDataFrame", data = meta),
              featureData = new("AnnotatedDataFrame", data = fmeta))
  
  # Creating a cello 
  Cello <- setClass("Cello",
                    slots = c(
                      name = "character", # The name of cvis
                      idx = "numeric", # The index of global cds object
                      proj = "list", # The projections as a list of data frame
                      pmeta = "data.frame", # The local meta data
                      notes = "character" # Other information to display to the user
                    )
  )
  cello1 <- new("Cello", name = cello.name, idx = 1:ncol(mtx)) 
  
  my_pca_proj <- seurat.obj@reductions$pca@cell.embeddings
  my_tsne_proj30 <- seurat.obj@reductions$tsne@cell.embeddings
  my_umap_proj30 <- seurat.obj@reductions$umap@cell.embeddings 
  
  cello1@proj <- list("PCA" = my_pca_proj, 
                      't-SNE30' = my_tsne_proj30,
                      "UMAP30" = my_umap_proj30) 
  
  if(extraDR){
    seurat.obj <- RunPCA(seurat.obj, npcs = 100, verbose = F,
                         assay = assay)
    seurat.obj <- regress_on_pca(seurat.obj, reg.var = paste0('nCount_ATAC'))
    
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:50, check_duplicates = FALSE,
                          assay = assay)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:50, verbose = F,
                          assay = assay)
    my_tsne_proj50 <- seurat.obj@reductions$tsne@cell.embeddings
    my_umap_proj50 <- seurat.obj@reductions$umap@cell.embeddings 
    
    
    
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:20, check_duplicates = FALSE,
                          assay = assay)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:20, verbose = F,
                          assay = assay)
    my_tsne_proj20 <- seurat.obj@reductions$tsne@cell.embeddings
    my_umap_proj20 <- seurat.obj@reductions$umap@cell.embeddings 
    
   
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:100, check_duplicates = FALSE,
                          assay = assay)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:100, verbose = F,
                          assay = assay)
    my_tsne_proj100 <- seurat.obj@reductions$tsne@cell.embeddings
    my_umap_proj100 <- seurat.obj@reductions$umap@cell.embeddings 
    
    cello1@proj <- list("PCA" = my_pca_proj, "t-SNE20" = my_tsne_proj20, 
                        "UMAP20" = my_umap_proj20, 't-SNE30' = my_tsne_proj30,
                        "UMAP30" = my_umap_proj30, 't-SNE50' = my_tsne_proj50,
                        "UMAP50" = my_umap_proj50, 't-SNE100' = my_tsne_proj100,
                        "UMAP100" = my_umap_proj100) 
    
    
  }
  
  clist <- list()
  clist[[cello.name]] <- cello1
  
  if(length(cluster_excluded) > 0){
    filteredCell = rownames(meta[!as.integer(as.character(meta$seurat_clusters)) %in% cluster_excluded, ])
    ids = which(colnames(seurat.obj) %in% filteredCell)
    cello2 <- new("Cello", name = cello.name, idx = ids) 
    
    
    cello2@proj <- list("PCA" = my_pca_proj[ids, ], 
                        't-SNE50' = my_tsne_proj50[ids, ],
                        "UMAP50" = my_umap_proj50[ids, ]) 
    if(extraDR){
      cello2@proj <- list("PCA" = my_pca_proj[ids, ], "t-SNE20" = my_tsne_proj20[ids, ], 
                          "UMAP20" = my_umap_proj20[ids, ], 't-SNE30' = my_tsne_proj30[ids, ],
                          "UMAP30" = my_umap_proj30[ids, ], 't-SNE50' = my_tsne_proj50[ids, ],
                          "UMAP50" = my_umap_proj50[ids, ], 't-SNE100' = my_tsne_proj100[ids, ],
                          "UMAP100" = my_umap_proj100[ids, ]) 
    }
    cello.name.clean = paste0(cello.name, '_clean')
    clist[[cello.name.clean]] <- cello2
  }
  
  
  return(list('eset' = eset, 'clist' = clist))
}





## prepare input4Cello using cisTopic model 
prepInput4Cello_cisTopic <- function(mtx, norm_mtx, cell_topic, bestModel = 100,
                                     extraMeta = NULL, cello.name = 'scATAC_cisTopic', 
                                     resl4clust = 0.3, 
                                     downSample = NULL){
  
  # create seurat to faciliate louvain clustering 
  seurat.obj <- CreateSeuratObject(mtx, assay = 'ATAC')
  seurat.obj@assays$ATAC@data <- norm_mtx
  
  kk = sapply(cell_topic, ncol)
  if(is.null(bestModel)) bestModel = max(kk)
  seurat.obj[['lda']] <- CreateDimReducObject(cell_topic[[which(kk == bestModel)]], key = 'Topic_',
                                              assay = DefaultAssay(seurat.obj))
  
  
  
  seurat.obj <- FindNeighbors(seurat.obj, reduction = 'lda', dims = 1:bestModel)
  seurat.obj <- FindClusters(seurat.obj, resolution = resl4clust)
  
  
  my_topic_proj <- cell_topic[[which(kk == bestModel)]]
  
  mdata = seurat.obj@meta.data
  fmeta <- data.frame(symbol = rownames(mtx)) 
  rownames(fmeta) <- fmeta$symbol
  
  # add extra meta data
  if(!is.null(extraMeta)){
    if(all(rownames(extraMeta) == rownames(mdata))) {
      shared.features = intersect(colnames(mdata), colnames(extraMeta))
      mdata[, shared.features] <- NULL
      mdata = cbind(mdata, extraMeta)
    }
  }
  
  # downsample
  ids = 1:ncol(seurat.obj)
  if(!is.null(downSample)){
    set.seed(2019)
    if(downSample <= ncol(seurat.obj)) 
      ids = sample(1:ncol(seurat.obj), downSample)
  }
  
  # Creating eset
  eset <- new("ExpressionSet",
              assayData = assayDataNew("environment", 
                                       exprs = mtx[, ids], 
                                       norm_exprs = norm_mtx[, ids]),
              phenoData =  new("AnnotatedDataFrame", data = mdata[ids, ]),
              featureData = new("AnnotatedDataFrame", data = fmeta))
  
  
  # Creating a cello 
  Cello <- setClass("Cello",
                    slots = c(
                      name = "character", # The name of cvis
                      idx = "numeric", # The index of global cds object
                      proj = "list", # The projections as a list of data frame
                      pmeta = "data.frame", # The local meta data
                      notes = "character" # Other information to display to the user
                    )
  )
  
  cello1 <- new("Cello", name = cello.name, idx = 1:length(ids)) 
  
  
  
  # cal tnse and umap for extra dims for cello
  for(dim0 in kk){
    tmp_topic = cell_topic[[which(kk == dim0)]]
    my_umap_proj0 <- uwot::umap(tmp_topic)
    my_tsne_proj0 <- Rtsne::Rtsne(tmp_topic, check_dup = F, pca = F)$Y
    rownames(my_tsne_proj0) = rownames(my_umap_proj0) <- rownames(cell_topic)
    colnames(my_umap_proj0) = c('UMAP_1', 'UMAP_2')
    colnames(my_tsne_proj0) = c('tSNE_1', 'tSNE_2')
    cello1@proj[[paste0('t-SEN', dim0)]] <- my_tsne_proj0[ids, ]
    cello1@proj[[paste0('UMAP', dim0)]] <- my_umap_proj0[ids, ]
    
  }
  
  clist <- list()
  clist[[cello.name]] <- cello1
  
  return(list('eset' = eset, 'clist' = clist))
}



## add extra meta data for each cell, with/without extra clustering
updateInput4Cello_cisTopic <- function(celloInputPath, cell_topic = NULL, 
                               extraMeta = NULL, clusterOnDim = NULL,
                               resolutions = c(0.2, 0.4, 0.6)){
  eset = readRDS(paste0(celloInputPath, '/eset.rds'))
  mdata = phenoData(eset)
  
  
  
  
  ## try different clustering
  if(!is.null(clusterOnDim) & !is.null(cell_topic)){
    kk = sapply(cell_topic, ncol)
    if(clusterOnDim %in% kk){
      sele_topic = cell_topic[[which(kk == clusterOnDim)]]
      mdata@data[[paste0('ATAC_snn_', clusterOnDim , '_res.', res0)]] <- 'something'
    }
      
  }
  
  if(!is.null(extraMeta)){
    if(all(rownames(extraMeta) == rownames(mdata))) {
      shared.features = intersect(colnames(mdata@data), colnames(extraMeta))
      mdata@data[, shared.features] <- NULL
      mdata@data = cbind(mdata@data, extraMeta)
    }
  }
  
  phenoData(eset) <- mdata
  
  
  saveRDS(eset, paste0(celloInputPath, '/eset.rds'))
}


## update mtx of cello objects, without touch other stuff
## need provide a mtx you want to add
updateRow4Cello_cisTopic <- function(celloInputPath, add_mtx, add_norm_mtx){
  eset = readRDS(paste0(celloInputPath, '/eset.rds'))
  mdata = phenoData(eset)
  mtx0 = eset@assayData$exprs
  norm_mtx0 = eset@assayData$norm_exprs
  
  add_mtx = add_mtx[, colnames(mtx0)]
  add_norm_mtx = add_norm_mtx[, colnames(mtx0)]
  
  mtx = rbind(mtx0, add_mtx)
  norm_mtx = rbind(norm_mtx0, add_norm_mtx)
  
  meta = mdata@data
  fmeta <- data.frame(symbol = rownames(mtx)) 
  rownames(fmeta) <- fmeta$symbol
  
  # preparing expression matrix
  eset <- new("ExpressionSet",
              assayData = assayDataNew("environment", 
                                       exprs = mtx, 
                                       norm_exprs = norm_mtx),
              phenoData =  new("AnnotatedDataFrame", data = meta),
              featureData = new("AnnotatedDataFrame", data = fmeta))
  
  
  phenoData(eset) <- mdata
  
  saveRDS(eset, paste0(celloInputPath, '/eset.rds'))
}



do_DA <- function(mtx_score, clusters, test = 'wilcox', fdr = 0.05, topn = 10){
  clusters$cluster = as.character(clusters$cluster)
  cls = unique(clusters$cluster)
  res = NULL
  features = rownames(mtx_score)
  for(cluster0 in cls){
    bc0 = clusters[cluster == cluster0]$barcode
    mtx1 = mtx_score[, colnames(mtx_score) %in% bc0]
    mtx2 = mtx_score[, !colnames(mtx_score) %in% bc0]
    mu1 = sapply(1:length(features), function(x) mean(mtx1[x, ]))
    mu2 = sapply(1:length(features), function(x) mean(mtx2[x, ]))
    
    
    pvs = sapply(1:length(features), function(x) wilcox.test(mtx1[x, ], mtx2[x, ], 
                                                             alternative = 'greater')$p.value )
    pvs.adj = p.adjust(pvs, method = 'fdr')
    res0 = data.table('feature' = features, 'cluster' = cluster0,
                      'mean1' = mu1, 'mean2' = mu2,
                      'pv' = pvs, 'pv_adjust' = pvs.adj)
    
    
    res0 = res0[order(pv_adjust), ]
    res0 = res0[pv_adjust <= fdr]
    
    if(nrow(res0) > topn) res0 = res0[1:topn, ]
    res = rbind(res, res0)
  }
  return(res)
}


#fg_genes: vector of forground genes
#bg_genes: background genes
#type: BP, CC, kegg
do_GO <- function(fg_genes, bg_genes = NULL, type = "BP", qCutoff = 0.05,
                  organism = c("mmu",  "hsa")) {
  if(organism =="mmu") {
    orgdb <- "org.Mm.eg.db"
    fromType = "SYMBOL"
    if(!require("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
    
  } else if(organism == "hsa") {
    orgdb <- "org.Hs.eg.db"
    if(!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
    fromType = "SYMBOL"
    
  }
  
  if(!is.null(bg_genes)) bg.df <- bitr(bg_genes, fromType = fromType,
                toType = c("SYMBOL", "ENTREZID"),
                OrgDb = orgdb)
  gene.df <- bitr(fg_genes, fromType = fromType,
                  toType = c("SYMBOL", "ENTREZID"),
                  OrgDb = orgdb)
  
  
  if(type == "kegg") {
    kegg_gene <- gene.df$ENTREZID
    if(!is.null(bg_genes)){
      kegg_bg <- bg.df$ENTREZID
      enrich_list <- enrichKEGG(
        gene          = kegg_gene,
        universe      = kegg_bg,
        organism      = organism,
        pAdjustMethod = "BH",
        qvalueCutoff  = qCutoff)
      
    }else{
      enrich_list <- enrichKEGG(
        gene          = kegg_gene,
        organism      = organism,
        pAdjustMethod = "BH",
        qvalueCutoff  = qCutoff)
      
    }
    
  }else {
    if(!is.null(bg_genes)){
      enrich_list <- enrichGO(gene        = gene.df$ENTREZID,
                              universe      = bg.df$ENTREZID,
                              OrgDb         = orgdb,
                              ont           = type,
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = qCutoff,
                              readable      = TRUE)
    }else{
      enrich_list <- enrichGO(gene        = gene.df$ENTREZID,
                              OrgDb         = orgdb,
                              ont           = type,
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = qCutoff,
                              readable      = TRUE)
    }
    
  }
  
  #return(enrich_list@result[enrich_list@result$qvalue <= qCutoff, ])
  return(enrich_list)
}

  
generate_gene_cisActivity <- function(gene_gtf, mtx, include_body = T,
                                      dist_to_tss = 2000){
  ## generating gene cis activity score
  ## input: gene gtf file, and atac matrix file
  
  ## 1. select gene up/down-stream regions (promoter with/without gene_body) ##
  
  gene_ann = fread(gene_gtf, sep = '\t')
  gene_ann = gene_ann[V3 == 'gene']
  gene_ann[, 'gene_name' := unlist(strsplit(V9, ';'))[3], by = V9]
  gene_ann[, 'gene_name' := gsub("\"", "", gene_name), by = gene_name]
  gene_ann[, 'gene_name' := unlist(strsplit(gene_name, ' '))[3], by = gene_name]
  names(gene_ann)[1] = 'chr'
  gene_ann = subset(gene_ann, select = c(chr, V4, V5, V7, gene_name))
  chrs = 1:22
  chrs = c(chrs, 'X', 'Y')
  gene_ann = gene_ann[chr %in% chrs]
  gene_ann = gene_ann[!duplicated(gene_name)]
  
  
  names(gene_ann)[2:4] = c('ss', 'ee', 'strand')
  if(!include_body){
    gene_ann[, 'start' := ifelse(strand == '+', ss - dist_to_tss, ee - dist_to_tss)]
    gene_ann[, 'end' := ifelse(strand == '+', ss + dist_to_tss, ee + dist_to_tss)]
  }else{
    gene_ann[, 'start' := ifelse(strand == '+', ss - dist_to_tss, ss)]
    gene_ann[, 'end' := ifelse(strand == '+', ee, ee + dist_to_tss)]
    
  }
  
  gene_ann[, 'chr' := paste0('chr', chr)]
  
  gene_ann = subset(gene_ann, select = c(chr, start, end, gene_name))
  
  
  ## 2. read mtx file ##
  
  rnames = rownames(mtx)
  chrs = sapply(rnames, function(x) unlist(strsplit(x, '-'))[1])
  starts = sapply(rnames, function(x) unlist(strsplit(x, '-'))[2])
  ends = sapply(rnames, function(x) unlist(strsplit(x, '-'))[3])
  
  peaks = data.table('chr' = chrs, 'start' = as.integer(starts), 
                     'end' = as.integer(ends))
  setkey(peaks, chr, start, end)
  peaks[, 'pname' := paste(chr, start, end, sep = '-')]
  over.ids = foverlaps(gene_ann, peaks, by.x = c('chr', 'start', 'end'),
                       by.y = c('chr', 'start', 'end'), which = T)
  over.ids[, 'gene_name' := gene_ann[xid, gene_name]]
  over.ids[, 'pname' := peaks[yid, pname]]
  over.ids = over.ids[complete.cases(over.ids)]
  
  
  smtx = sparseMatrix(i = over.ids$xid, j = over.ids$yid,
                      dimnames = list(gene_ann$gene_name[1:max(over.ids$xid)],
                                      peaks$pname[1:max(over.ids$yid)]))
  
  mtx = mtx[rownames(mtx) %in% colnames(smtx), ]
  smtx = smtx[, rownames(mtx)]
  
  activity.matrix = smtx %*% mtx
  rs = Matrix::rowSums(activity.matrix)
  activity.matrix = activity.matrix[rs > 10, ]
  
  return(activity.matrix)
  
}


## informative genes selection using gini index from qin
## if gini_cut_qt is an integer, it will output top 1:gini_cut_qt genes as include genes
ifg_select <- function(data, cluster, cluster_min_cell_num = 100, 
                       min_cluster_expr_fraction = .1, gini_cut_qt = .75, 
                       do_go = F, filePath = NULL, 
                       fileName = "ifg_select_",
                       orgdb = "org.Hs.eg.db", gene_id_type = "SYMBOL",
                       gene_background = NULL) {
  
  if(ncol(data) != length(cluster)) stop("Cell number do not match cluster length.")
  
  use_clus <- names(which(table(cluster) >= cluster_min_cell_num))
  
  # For each cluster, compute gene expressed fraction
  
  expr_clus_frac <- sapply(use_clus, function(x) {
    
    cur_data <- data[, cluster == x]
    
    Matrix::rowMeans(cur_data > 0)
    
  })
  
  # Compute gini coefficient 
  
  # Require a gene to be expressed in at least one cluster with at least .1 expressed fraction to be considered for downstream uniform gene selection
  
  use_g <- rownames(expr_clus_frac)[rowSums(expr_clus_frac >= min_cluster_expr_fraction) > 0] # 11813
  
  message(paste0("Selecting informative features from ", length(use_g), " robustly detected features."))
  
  expr_clus_frac <- expr_clus_frac[rownames(expr_clus_frac) %in% use_g,]
  
  gene_clus_gini <- apply(expr_clus_frac, 1, gini)
  
  
  
  #pdf(paste0(filePath, fileName, "gene_clus_gini_hist.pdf"))
  
  #hist(gene_clus_gini, breaks = 100)
  
  #dev.off()
  
  
  
  gene_clus_gini = sort(gene_clus_gini, decreasing = T)
  if(gini_cut_qt > 1){
    topn = min(length(use_g), gini_cut_qt)
    exclude_g <- names(gene_clus_gini[-(1:topn)])
    
    include_g <- names(gene_clus_gini[(1:topn)])
  }else{
    gini_cut <- quantile(gene_clus_gini, gini_cut_qt) 
    
    message(paste0("Cut at gini quantile ", gini_cut_qt, " with value ", gini_cut))
    
    exclude_g <- names(gene_clus_gini)[gene_clus_gini < gini_cut]
    
    include_g <- names(gene_clus_gini)[gene_clus_gini >= gini_cut]
    
  }
  
  #write.csv(data.frame(less_specific_feature = exclude_g), paste0(filePath, fileName, "less_specific_feature_list.csv"))
  
  #write.csv(data.frame(specific_feature = include_g), paste0(filePath, fileName, "specific_feature_list.csv"))
  
  
  
  message(paste0("Found ", length(exclude_g), " less specific features."))
  
  message(paste0("Returning ", length(include_g), " specific features."))
  
  
  
  if(do_go) {
    
    if(!length(gene_background)) {
      
      gene_background <- use_g
      
    }
    
    exclude_g.df <- bitr(exclude_g, fromType = gene_id_type,
                         
                         toType = c("ENTREZID"),
                         
                         OrgDb = orgdb)
    
    bg.df <- bitr(gene_background, fromType = gene_id_type,
                  
                  toType = c("ENTREZID"),
                  
                  OrgDb = orgdb)
    
    exclude_g_go<- enrichGO(gene        = exclude_g.df$ENTREZID,
                            
                            universe      = bg.df$ENTREZID,
                            
                            OrgDb         = orgdb,
                            
                            ont           = "BP",
                            
                            pAdjustMethod = "BH",
                            
                            pvalueCutoff  = 0.01,
                            
                            qvalueCutoff  = 0.05,
                            
                            readable      = TRUE)
    
    write.csv(exclude_g_go, paste0(filePath, fileName, "exlude_feature_go.csv"))
    
    
    
    
    
    include_g.df <- bitr(include_g, fromType = gene_id_type,
                         
                         toType = c("ENTREZID"),
                         
                         OrgDb = orgdb)
    
    bg.df <- bitr(gene_background, fromType = gene_id_type,
                  
                  toType = c("ENTREZID"),
                  
                  OrgDb = orgdb)
    
    include_g_go <- enrichGO(gene        = include_g.df$ENTREZID,
                            
                            universe      = bg.df$ENTREZID,
                            
                            OrgDb         = orgdb,
                            
                            ont           = "BP",
                            
                            pAdjustMethod = "BH",
                            
                            pvalueCutoff  = 0.01,
                            
                            qvalueCutoff  = 0.05,
                            
                            readable      = TRUE)
    
    write.csv(include_g_go, paste0(filePath, fileName, "include_feature_go.csv"))
    
  }
  
  
  
  return(list('include_g' = include_g, 'exclude_g' = exclude_g))
  
}


## re-select genes/features by F-stat
fstat_select <- function(data, cluster, cluster_min_cell_num = 100, 
                       min_cluster_expr_fraction = .1, cut_qt = 0.75, 
                       do_go = F, filePath = NULL, 
                       fileName = "fstat_select_",
                       orgdb = "org.Hs.eg.db", gene_id_type = "SYMBOL",
                       gene_background = NULL,
                       max_cell_per_cl = 200) {
  
  if(ncol(data) != length(cluster)) stop("Cell number do not match cluster length.")
  nn = names(cluster)
  cluster = as.character(cluster)
  names(cluster) = nn
  
  use_clus <- names(which(table(cluster) >= cluster_min_cell_num))
  
  # For each cluster, compute gene expressed fraction
  
  expr_clus_frac <- sapply(use_clus, function(x) {
    
    cur_data <- data[, cluster == x]
    
    Matrix::rowMeans(cur_data > 0)
    
  })
  
  # Compute F-stats
  
  # Require a gene to be expressed in at least one cluster with at least .1 expressed fraction to be considered for downstream uniform gene selection
  
  use_g <- rownames(expr_clus_frac)[rowSums(expr_clus_frac >= min_cluster_expr_fraction) > 0] 
  
  message(paste0("Selecting informative features from ", length(use_g), " robustly detected features."))
  
  expr_clus_frac <- expr_clus_frac[rownames(expr_clus_frac) %in% use_g,]
  
  data = data[rownames(data) %in% use_g, ]
  data = log1p(data) 
  
  cluster = cluster[cluster %in% use_clus]
  data = data[, names(cluster)]
  
  ##downsample each cluster to have at most 200 cells
  cluster0 = NULL
  set.seed(2019)
  for(cl in use_clus){
    tmp = cluster[cluster == cl]
    if(length(tmp) > max_cell_per_cl){
      cluster0 = c(cluster0, tmp[sample(1:length(tmp), max_cell_per_cl)])
    }else{
      cluster0 = c(cluster0, tmp)
    }
  }
  
  data = data[, names(cluster0)]
  data = (data > 0) * 1
  gene_clus_fstat <- apply(data, 1, function(x) anova(lm(x ~ cluster0))$`F value`[1])
  #gene_clus_fstat = rep(0, nrow(data))
  #for(i in 1:nrow(data)){
  #  gene_clus_fstat[i] = anova(lm(data[i, ] ~ cluster))$`F value`[1]
  #}
  
  #gene_clus_fstat = pmin(gene_clus_fstat * length(gene_clus_fstat), 1)
  
  pdf(paste0(filePath, fileName, "gene_clus_fstat_pv_hist.pdf"))
  
  hist(gene_clus_fstat, breaks = 100)
  
  dev.off()
  
  cut_thr = quantile(gene_clus_fstat, cut_qt)
  exclude_g <- names(gene_clus_fstat)[gene_clus_fstat < cut_thr]
  
  include_g <- names(gene_clus_fstat)[gene_clus_fstat >= cut_thr]
  
  #write.csv(data.frame(less_specific_feature = exclude_g), paste0(filePath, fileName, "less_specific_feature_list.csv"))
  
  #write.csv(data.frame(specific_feature = include_g), paste0(filePath, fileName, "specific_feature_list.csv"))
  
  
  
  message(paste0("Found ", length(exclude_g), " less specific features."))
  
  message(paste0("Returning ", length(include_g), " specific features."))
  
  
  
  if(do_go) {
    
    if(!length(gene_background)) {
      
      gene_background <- use_g
      
    }
    
    exclude_g.df <- bitr(exclude_g, fromType = gene_id_type,
                         
                         toType = c("ENTREZID"),
                         
                         OrgDb = orgdb)
    
    bg.df <- bitr(gene_background, fromType = gene_id_type,
                  
                  toType = c("ENTREZID"),
                  
                  OrgDb = orgdb)
    
    exclude_g_go<- enrichGO(gene        = exclude_g.df$ENTREZID,
                            
                            universe      = bg.df$ENTREZID,
                            
                            OrgDb         = orgdb,
                            
                            ont           = "BP",
                            
                            pAdjustMethod = "BH",
                            
                            pvalueCutoff  = 0.01,
                            
                            qvalueCutoff  = 0.05,
                            
                            readable      = TRUE)
    
    write.csv(exclude_g_go, paste0(filePath, fileName, "exlude_feature_go.csv"))
    
    
    
    
    
    include_g.df <- bitr(include_g, fromType = gene_id_type,
                         
                         toType = c("ENTREZID"),
                         
                         OrgDb = orgdb)
    
    bg.df <- bitr(gene_background, fromType = gene_id_type,
                  
                  toType = c("ENTREZID"),
                  
                  OrgDb = orgdb)
    
    include_g_go <- enrichGO(gene        = include_g.df$ENTREZID,
                             
                             universe      = bg.df$ENTREZID,
                             
                             OrgDb         = orgdb,
                             
                             ont           = "BP",
                             
                             pAdjustMethod = "BH",
                             
                             pvalueCutoff  = 0.01,
                             
                             qvalueCutoff  = 0.05,
                             
                             readable      = TRUE)
    
    write.csv(include_g_go, paste0(filePath, fileName, "include_feature_go.csv"))
    
  }
  
  
  
  return(list('include_g' = include_g, 'exclude_g' = exclude_g))
  
}


## undate seurat given variable features
updateSeurat_gvFeatures <- function(seurat.obj, gvFeatures, reg.var = 'nCount_RNA',
                                  doScale = T, doCenter = T, npc = 50,
                                  clusterOn = 'pca', resolution = 0.5){
  
  
  VariableFeatures(seurat.obj) <- gvFeatures
 
  seurat.obj <- ScaleData(object = seurat.obj, 
                            features = VariableFeatures(seurat.obj), 
                            vars.to.regress = reg.var, do.scale = doScale, do.center = doCenter)
    
    
 
  seurat.obj <- RunPCA(object = seurat.obj, 
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  seurat.obj <- RunTSNE(object = seurat.obj, reduction = 'pca', dims = 1:npc,
                        check_duplicates = FALSE)
  seurat.obj <- RunUMAP(object = seurat.obj, 
                       reduction = 'pca',
                       verbose = FALSE, dims = 1:npc)
  seurat.obj <- FindNeighbors(object = seurat.obj, reduction = clusterOn)
  seurat.obj <- FindClusters(object = seurat.obj,  resolution = resolution)
  
  return(seurat.obj)
}
updateSeurat_gvFeatures = cmpfun(updateSeurat_gvFeatures)




#integrate a list of seurat objects
#supporse each element of the list was normalized and/or scaled
integrateSeuratList <- function(seurat.list, npc = 50, reg.var = NULL,
                               reduction = 'pca', anchor.features = 6000){
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                           anchor.features = anchor.features)
  rm(seurat.list)
  
  seurat.integrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:npc)
  rm(seurat.anchors)
  
  DefaultAssay(seurat.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE, vars.to.regress = reg.var)
  seurat.integrated <- RunPCA(seurat.integrated, npcs = npc, verbose = FALSE)
  #seurat.integrated <- RunTSNE(seurat.integrated, reduction = reduction, dims = 1:npc,
  #                             check_duplicates = FALSE)
  #seurat.integrated <- RunUMAP(seurat.integrated, reduction = reduction, dims = 1:npc)
 
  return(seurat.integrated)
  
}

#integrate a list of seurat objects
integrateSeuratList_withReference <- function(seurat.list, npc = 50, reg.var = NULL,
                                reduction = 'pca', anchor.features = 6000,
                                ref.name = 'sample1'){
  ref.dataset <- which(names(seurat.list) == ref.name)
  
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                           anchor.features = anchor.features,
                                           reference = ref.dataset)
  rm(seurat.list)
  
  seurat.integrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:npc)
  rm(seurat.anchors)
  
  DefaultAssay(seurat.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE, vars.to.regress = reg.var)
  seurat.integrated <- RunPCA(seurat.integrated, npcs = npc, verbose = FALSE)
  #seurat.integrated <- RunTSNE(seurat.integrated, reduction = reduction, dims = 1:npc,
  #                             check_duplicates = FALSE, 
  #                             tsne.method = "FIt-SNE", nthreads = 2, max_iter = 2000)
  #seurat.integrated <- RunUMAP(seurat.integrated, reduction = reduction, dims = 1:npc)
  
  return(seurat.integrated)
  
}



#integrate a list of seurat objects by SCT normalization
integrateSeuratList_SCT <- function(seurat.list, npc = 50, reg.var = NULL,
                                reduction = 'pca', anchor.features = 3000,
                                resolution = 0.6){
  anchor.features = SelectIntegrationFeatures(seurat.list, 
                                            nfeatures = anchor.features)
  seurat.list <- PrepSCTIntegration(object.list = seurat.list, 
                                           verbose = F,
                                           anchor.features = anchor.features)
  seurat.anchors = FindIntegrationAnchors(object.list = seurat.list, 
                                          normalization.method = "SCT", 
                                          anchor.features = anchor.features, 
                                          verbose = FALSE)
  rm(seurat.list)
  
  seurat.integrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:npc,
                                     normalization.method = "SCT", verbose = F)
  rm(seurat.anchors)
  
  DefaultAssay(seurat.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  seurat.integrated <- RunPCA(seurat.integrated, npcs = npc, verbose = FALSE)
  seurat.integrated <- RunTSNE(seurat.integrated, reduction = reduction, dims = 1:npc)
  seurat.integrated <- RunUMAP(seurat.integrated, reduction = reduction, dims = 1:npc)
  seurat.integrated <- FindNeighbors(seurat.integrated, reduction = reduction, 
                                     dims = 1:npc)
  seurat.integrated <- FindClusters(seurat.integrated, resolution = resolution)
  
  return(seurat.integrated)
  
}


## define ctype signature score by marker genes
ctype_assign_chisq <- function(mtx, ctype_markers, 
                               upper_thr = 0.1, lower_thr = 0.5,
                               negative_marker_w = 0.3){
  ctypes = unique(ctype_markers$ctype)
  
  mtx = as.matrix(mtx[rownames(mtx) %in% unique(ctype_markers$gene), ])
  score.mtx = chisq.test(mtx)$residuals
  score.mtx[is.na(score.mtx)] = 0
  
  score.mtx = as.matrix(score.mtx[rownames(score.mtx) %in% unique(ctype_markers$gene), ])
  score.mtx.neg = score.mtx.pos = score.mtx
  for(i in 1:nrow(score.mtx)){
    tmp = score.mtx[i,]
    uthr = quantile(tmp[tmp > 0], upper_thr)
    tmp_pos = ifelse(tmp >= uthr, 1, 0)
    dthr = quantile(-tmp[tmp < 0], lower_thr)
    tmp_neg = ifelse(tmp <= -dthr, 1, 0)
    
    score.mtx.pos[i, ] = tmp_pos
    score.mtx.neg[i, ] = tmp_neg
    
    
  }
  
  # calculate support freq for each marker gene given a ctype
  supp_freq = list()
  score_cell = list()
  for(ctype0 in ctypes){
    genes0_freq = ctype_markers[ctype == ctype0, ]
    genes0_freq = genes0_freq[order(-N), ]
    if(length(which(rownames(score.mtx) %in% genes0_freq$gene)) == 0) next
    genes0_freq = genes0_freq[gene %in% rownames(score.mtx), ]
    
    # different weight for neg and postive markers
    genes0_freq[, freq := ifelse(direction == -1, 
                                 freq * negative_marker_w, freq *(1-negative_marker_w))]
    genes0_freq$freq = genes0_freq$freq/sum(genes0_freq$freq)
    genes0_freq.pos = genes0_freq[direction == 1]
    genes0_freq.neg = genes0_freq[direction == -1]
    
    
    score.mtx0.pos = score.mtx.pos[rownames(score.mtx.pos) %in% genes0_freq.pos$gene, ]
    
    if(nrow(genes0_freq.pos) == 1){
      score_cell[[ctype0]] = score.mtx0.pos * genes0_freq.pos$freq 
    }else{
      score_cell[[ctype0]] = t(genes0_freq.pos$freq) %*% score.mtx0.pos
    }
    
    if(any(genes0_freq$direction == -1)){
      score.mtx0.neg = score.mtx.neg[rownames(score.mtx.pos) %in% genes0_freq.neg$gene, ]
      
      if(nrow(genes0_freq.neg) == 1){
        score_cell[[ctype0]] = score_cell[[ctype0]] + score.mtx0.neg * genes0_freq.neg$freq 
      }else{
        score_cell[[ctype0]] = score_cell[[ctype0]] + t(genes0_freq.neg$freq) %*% score.mtx0.neg
      }
    }
    
    
    
    
  }
  
  res = do.call('rbind', score_cell)
  rownames(res) = names(score_cell)
  return(res)
}


## define ctype signature score by marker genes
ctype_assign_chisq_update <- function(mtx, ctype_markers, 
                               upper_thr = 0.1, lower_thr = 0.1,
                               negative_marker_w = 0.5){
  ctypes = unique(ctype_markers$ctype)
  
  mtx = as.matrix(mtx[rownames(mtx) %in% unique(ctype_markers$gene), ])
  score.mtx = chisq.test(mtx)$residuals
  score.mtx[is.na(score.mtx)] = 0
  
  score.mtx = as.matrix(score.mtx[rownames(score.mtx) %in% unique(ctype_markers$gene), ])
  score.mtx.neg = score.mtx.pos = score.mtx
  for(i in 1:nrow(score.mtx)){
    tmp = score.mtx[i,]
    uthr = quantile(tmp[tmp > 0], upper_thr)
    tmp_pos = ifelse(tmp >= uthr, 1, 0)
    dthr = quantile(tmp[tmp > 0], lower_thr)
    tmp_neg = ifelse(tmp <= dthr, 1, 0)
    
    score.mtx.pos[i, ] = tmp_pos
    score.mtx.neg[i, ] = tmp_neg
    
    
  }
  
  # calculate support freq for each marker gene given a ctype
  supp_freq = list()
  score_cell = list()
  for(ctype0 in ctypes){
    genes0_freq = ctype_markers[ctype == ctype0, ]
    genes0_freq = genes0_freq[order(-N), ]
    if(length(which(rownames(score.mtx) %in% genes0_freq$gene)) == 0) next
    genes0_freq = genes0_freq[gene %in% rownames(score.mtx), ]
    
    # different weight for neg and postive markers
    genes0_freq[, freq := ifelse(direction == -1, 
                                 freq * negative_marker_w, freq *(1-negative_marker_w))]
    genes0_freq$freq = genes0_freq$freq/sum(genes0_freq$freq)
    genes0_freq.pos = genes0_freq[direction == 1]
    genes0_freq.neg = genes0_freq[direction == -1]
    
    
    score.mtx0.pos = score.mtx.pos[rownames(score.mtx.pos) %in% genes0_freq.pos$gene, ]
    
    if(nrow(genes0_freq.pos) == 1){
      score_cell[[ctype0]] = score.mtx0.pos * genes0_freq.pos$freq 
    }else{
      score_cell[[ctype0]] = t(genes0_freq.pos$freq) %*% score.mtx0.pos
    }
    
    if(any(genes0_freq$direction == -1)){
      score.mtx0.neg = score.mtx.neg[rownames(score.mtx.pos) %in% genes0_freq.neg$gene, ]
      
      if(nrow(genes0_freq.neg) == 1){
        score_cell[[ctype0]] = score_cell[[ctype0]] + score.mtx0.neg * genes0_freq.neg$freq 
      }else{
        score_cell[[ctype0]] = score_cell[[ctype0]] + t(genes0_freq.neg$freq) %*% score.mtx0.neg
      }
    }
    
    
    
    
  }
  
  res = do.call('rbind', score_cell)
  rownames(res) = names(score_cell)
  return(res)
}



signatureScore_zscore <- function(seurat.obj, pfeatures, nfeatures = NULL,
                                  vars.to.regress = NULL,
                                  score.name = 'quicent'){
  
  pfeatures = pfeatures[pfeatures %in% rownames(seurat.obj)]
  if(!is.null(nfeatures)) nfeatures = nfeatures[nfeatures %in% rownames(seurat.obj)]
  
  seurat.obj <- ScaleData(seurat.obj, features = c(pfeatures, nfeatures),
                          do.scale = T, do.center = T, 
                          vars.to.regress = vars.to.regress)
  
  mtx = GetAssayData(seurat.obj, slot = 'scale.data')
  
  pscore = mtx[pfeatures, ]
  
  if(length(nfeatures) > 0) nscore = mtx[nfeatures, ]
  
  if(length(pfeatures) > 1) pscore = Matrix::colSums(pscore)
  if(length(nfeatures) > 1)  nscore = Matrix::colSums(nscore)
  if(length(nfeatures) > 0){
    scores = (pscore - nscore)/length(c(pfeatures, nfeatures))
  }else{
    scores = pscore/length(pfeatures)
  }
  seurat.obj <- AddMetaData(seurat.obj, metadata = scores, 
                            col.name = score.name)
  return(seurat.obj)
}



## codes may not used any more ####
# rebin data in unit of tad, adjusting tad size

rebin_matrix2tad <- function(mtx, tads){
  # mtx: matrix wiht rownames as chr-start-end and
  # colnames as cell names
  tads[, 'size' := floor(end/1000 - start/1000)]
  setkey(tads, id)
  
  rnames = rownames(mtx)
  mtx_chr = sapply(rnames, function(x) unlist(strsplit(x, '-'))[1])
  chrs = unique(mtx_chr)
  starts = as.numeric(sapply(rnames, function(x) unlist(strsplit(x, '-'))[2]))
  ends = as.numeric(sapply(rnames, function(x) unlist(strsplit(x, '-'))[3]))
  
  rebin_mat = NULL
  for(chr0 in chrs){
    mtx0 = mtx[mtx_chr == chr0, ]
    mtx0 = as.data.table(mtx0)
    mtx0$start = starts[mtx_chr == chr0] 
    mtx0$end = ends[mtx_chr == chr0] 
    mtx0[, 'midP' := start/2 + end/2]
    tads0 = tads[chr == chr0]
    mtx0[, 'tad_id' := which(tads0$start <= midP & tads0$end >= midP), by = midP]
    mtx0[, tad_id := tads0$id[tad_id]]
    mtx0[, c('start', 'end', 'midP') := NULL]
    rebin_mat = rbind(rebin_mat, mtx0)
  }
  
  rebin_mat = data.table(rebin_mat)
  setkey(rebin_mat, tad_id)
  
  new_mat = rebin_mat[, lapply(.SD, sum), by = tad_id]
  new_mat = new_mat[complete.cases(new_mat)]
  
  feature.names = new_mat$tad_id
  new_mat[, 'size' := tads[J(new_mat$tad_id)]$size]
  
  
  new_mat = new_mat[, lapply(.SD, function(x) x/size * 1000),
                    .SDcols = !c('size', 'tad_id')]
  
  new_mat = as.matrix(new_mat)
  rownames(new_mat) = feature.names
  
  return(new_mat)
}
rebin_matrix2tad = cmpfun(rebin_matrix2tad)




# in tad
doBasicSeurat_tad <- function(tads, mtx, npc = 50){
  
  tads = full_seg_tads(tads)
  tads[, 'id' := paste(chr, start, end, sep = '-')]
  rebinned_mtx = rebin_matrix2tad(mtx, tads)
  rebinned_mtx = log2(1 + rebinned_mtx)
  seurat.obj = CreateSeuratObject(rebinned_mtx, project = 'scATAC_tad', assay = 'ATAC',
                                  names.delim = '-')
  cell.names = colnames(rebinned_mtx)
  
  #seurat.obj <- NormalizeData(seurat.obj, normalization.method = 'LogNormalize',
  #                            scale.factor = 1e4)
  
  seurat.obj <- FindVariableFeatures(object = seurat.obj, 
                                     selection.method = 'dispersion', 
                                     nfeatures = nrow(rebinned_mtx))
  seurat.obj <- ScaleData(object = seurat.obj, 
                          features = VariableFeatures(seurat.obj), 
                          vars.to.regress = c('nCount_ATAC'), do.scale = T, do.center = T)
  
  
  seurat.obj <- RunPCA(object = seurat.obj, 
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  
  return(seurat.obj)
}
doBasicSeurat_tad = cmpfun(doBasicSeurat_tad)

# in bin
doBasicSeurat_bin <- function(mtx, resl = 500 * 1000, npc = 50){
  
  rebinned_mtx = rebin_matrix2Bin(mtx, resl)
  rebinned_mtx = log2(1 + rebinned_mtx)
  seurat.obj = CreateSeuratObject(rebinned_mtx, project = 'scATAC_tad', assay = 'ATAC',
                                  names.delim = '-')
  cell.names = colnames(rebinned_mtx)
  
  seurat.obj <- FindVariableFeatures(object = seurat.obj, 
                                     selection.method = 'dispersion', 
                                     nfeatures = nrow(rebinned_mtx))
  seurat.obj <- ScaleData(object = seurat.obj, 
                          features = VariableFeatures(seurat.obj), 
                          vars.to.regress = c('nCount_ATAC'), do.scale = T, do.center = T)
  
  
  seurat.obj <- RunPCA(object = seurat.obj, 
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  
  return(seurat.obj)
}
doBasicSeurat_bin = cmpfun(doBasicSeurat_bin)


# in given ctcf peaks
doBasicSeurat_ctcf <- function(ctcf_peaks, mtx, npc = 50){
  
  ptads = segByPeak(ctcf_peaks)
  ptads[, 'id' := paste(chr, start, end, sep = '-')]
  ptads[, 'size' := floor((end - start)/1000)]
  ptads = ptads[size > 2 & size < 3000]
  rebinned_mtx = rebin_matrix2tad(mtx, ptads)
  rebinned_mtx = log2(1 + rebinned_mtx)
  seurat.obj = CreateSeuratObject(rebinned_mtx, project = 'scATAC_tad', assay = 'ATAC',
                                  names.delim = '-')
  cell.names = colnames(rebinned_mtx)
  
  #seurat.obj <- NormalizeData(seurat.obj, normalization.method = 'LogNormalize',
  #                            scale.factor = 1e4)
  
  seurat.obj <- FindVariableFeatures(object = seurat.obj, 
                                     selection.method = 'dispersion', 
                                     nfeatures = nrow(rebinned_mtx) * 0.5)
  seurat.obj <- ScaleData(object = seurat.obj, 
                          features = VariableFeatures(seurat.obj), 
                          vars.to.regress = c('nCount_ATAC'), do.scale = T, do.center = T)
  
  
  seurat.obj <- RunPCA(object = seurat.obj, 
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  
  return(seurat.obj)
}
doBasicSeurat_ctcf = cmpfun(doBasicSeurat_ctcf)


prepInput4Cello_from_liger <- function(liger.obj, cello.name = 'LIGER_OUT'){
  
  ## get union featured matrix
  
  len = length(liger.obj@raw.data)
  mtx = t(liger.obj@raw.data[[1]])
  norm.mtx = t(liger.obj@norm.data[[1]])
  if(len > 1){
    for(i in 2:len){
      bb = t(liger.obj@raw.data[[i]])
      cc = t(liger.obj@norm.data[[i]])
      mtx = Matrix.utils::rBind.fill(mtx, bb)
      norm.mtx = Matrix.utils::rBind.fill(norm.mtx, bb)
    }
  }
  
  mtx = t(mtx)
  norm.mtx = t(norm.mtx)
  
  
  meta = liger.obj@cell.data
  meta = cbind(meta, 'cluster' = liger.obj@clusters)
  
  fmeta <- data.frame(symbol = rownames(mtx)) 
  rownames(fmeta) <- fmeta$symbol
  
  eset <- new("ExpressionSet",
              assayData = assayDataNew("environment", 
                                       exprs = mtx, 
                                       norm_exprs = norm.mtx),
              phenoData =  new("AnnotatedDataFrame", data = meta),
              featureData = new("AnnotatedDataFrame", data = fmeta))
  
  # Creating a cello 
  Cello <- setClass("Cello",
                    slots = c(
                      name = "character", # The name of cvis
                      idx = "numeric", # The index of global cds object
                      proj = "list", # The projections as a list of data frame
                      pmeta = "data.frame", # The local meta data
                      notes = "character" # Other information to display to the user
                    )
  )
  cello1 <- new("Cello", name = cello.name, idx = 1:ncol(mtx)) 
  
  
  my_umap_proj <- liger.obj@tsne.coords 
  colnames(my_umap_proj) = c('UMAP_1', 'UMAP_2')
  cello1@proj <- list("UMAP" = my_umap_proj) 
  
  
  clist <- list()
  clist[[cello.name]] <- cello1
  
  
  return(list('eset' = eset, 'clist' = clist))
}


prepInput4Cello_gene_cis_from_liger <- function(liger.obj, cello.name = 'LIGER_OUT', gtf_file){
  
  ## get union featured matrix
  
  len = length(liger.obj@raw.data)
  mtx = t(liger.obj@raw.data[[1]])
  #norm.mtx = t(liger.obj@norm.data[[1]])
  if(len > 1){
    for(i in 2:len){
      bb = t(liger.obj@raw.data[[i]])
      #cc = t(liger.obj@norm.data[[i]])
      mtx = Matrix.utils::rBind.fill(mtx, bb, fill = 0)
      #norm.mtx = Matrix.utils::rBind.fill(norm.mtx, bb)
    }
  }
  
  mtx = t(mtx)
  #norm.mtx = t(norm.mtx)
  
  
  ## get gene cis score
  activity.mtx = generate_gene_cisActivity(gtf_file, mtx)
  rm(mtx)
  
  
  meta = liger.obj@cell.data
  meta = cbind(meta, 'cluster' = liger.obj@clusters)
  
  fmeta <- data.frame(symbol = rownames(activity.mtx)) 
  rownames(fmeta) <- fmeta$symbol
  
  eset <- new("ExpressionSet",
              assayData = assayDataNew("environment", 
                                       exprs = activity.mtx, 
                                       norm_exprs = log2(activity.mtx+1)),
              phenoData =  new("AnnotatedDataFrame", data = meta),
              featureData = new("AnnotatedDataFrame", data = fmeta))
  
  # Creating a cello 
  Cello <- setClass("Cello",
                    slots = c(
                      name = "character", # The name of cvis
                      idx = "numeric", # The index of global cds object
                      proj = "list", # The projections as a list of data frame
                      pmeta = "data.frame", # The local meta data
                      notes = "character" # Other information to display to the user
                    )
  )
  cello1 <- new("Cello", name = cello.name, idx = 1:ncol(activity.mtx)) 
  
  
  my_umap_proj <- liger.obj@tsne.coords 
  colnames(my_umap_proj) = c('UMAP_1', 'UMAP_2')
  cello1@proj <- list("UMAP" = my_umap_proj) 
  
  
  clist <- list()
  clist[[cello.name]] <- cello1
  
  
  return(list('eset' = eset, 'clist' = clist))
}

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, myColors,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) + scale_fill_manual(values = myColors)
  return(p)
}
## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}
## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, myColors,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, myColors = myColors,
                                                              pt.size = pt.size, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}



