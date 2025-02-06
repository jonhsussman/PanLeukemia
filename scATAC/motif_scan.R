library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(GenomicRanges)

## check motif overlapping with any peak -- AML ####

seurat.atac.AML <- readRDS(file = 'PanLeukemia/scATAC/Motif/AML_scATAC_downsampleBySeuratProjection_Signca_RecalledPeaks_BlacklistRemoved_Oct12_2022_withChromVAR.rds')

DefaultAssay(seurat.atac.AML) <- 'ATAC'

peaks = rownames(seurat.atac.AML)
peaks = sapply(peaks, function(x) unlist(strsplit(x, ','))[1])
names(peaks) = NULL

pks = tidyr::separate(data = data.table('peak_name' = peaks),
                      col = 'peak_name', into = c('chr', 'start', 'end'),
                      remove = F)
pks$start = as.integer(pks$start)
pks$end = as.integer(pks$end)
setkey(pks, chr, start)
# Make a set of peaks
peaks <- GenomicRanges::GRanges(seqnames = pks$chr,
                 ranges = IRanges::IRanges(start = pks$start,
                                  end = pks$end))

motif_ix <- matchMotifs(human_pwms_v2, peaks, 
                        genome = BSgenome.Hsapiens.UCSC.hg38,
                        p.cutoff = 1e-3)
motif_ix_str <- matchMotifs(human_pwms_v2, peaks, 
                        genome = BSgenome.Hsapiens.UCSC.hg38,
                        p.cutoff = 5e-5)

if(F){
        motif_name = motif_ix@colData$name
        pk_inf = motif_ix@rowRanges
        motif_pk_match <- motif_ix@assays@data$motifMatches
        rownames(motif_pk_match) = paste0(pk_inf@seqnames, '-', pk_inf@ranges)
        names(motif_name) = NULL
        colnames(motif_pk_match) = motif_name      
}

saveRDS(motif_ix, file = 'TRN/Motifs/AML_motif_match_pvE-03.rds')
saveRDS(motif_ix_str, file = 'TRN/Motifs/AML_motif_match_pv5E-05.rds')


## check motif overlapping with any peak -- B-ALL ####

seurat.atac.BALL <- readRDS(file = 'PanLeukemia/scATAC/Motif/BALL_scATAC_downsampleBySeuratProjection_Signca_RecalledPeaks_BlacklistRemoved_Oct12_2022_withChromVAR.rds')

DefaultAssay(seurat.atac.BALL) <- 'ATAC'

peaks = rownames(seurat.atac.BALL)
peaks = sapply(peaks, function(x) unlist(strsplit(x, ','))[1])
names(peaks) = NULL

pks = tidyr::separate(data = data.table('peak_name' = peaks),
                      col = 'peak_name', into = c('chr', 'start', 'end'),
                      remove = F)
pks$start = as.integer(pks$start)
pks$end = as.integer(pks$end)
setkey(pks, chr, start)
# Make a set of peaks
peaks <- GenomicRanges::GRanges(seqnames = pks$chr,
                                ranges = IRanges::IRanges(start = pks$start,
                                                          end = pks$end))

motif_ix <- matchMotifs(human_pwms_v2, peaks, 
                        genome = BSgenome.Hsapiens.UCSC.hg38,
                        p.cutoff = 1e-3)
motif_ix_str <- matchMotifs(human_pwms_v2, peaks, 
                            genome = BSgenome.Hsapiens.UCSC.hg38,
                            p.cutoff = 5e-5)

if(F){
  motif_name = motif_ix@colData$name
  pk_inf = motif_ix@rowRanges
  motif_pk_match <- motif_ix@assays@data$motifMatches
  rownames(motif_pk_match) = paste0(pk_inf@seqnames, '-', pk_inf@ranges)
  names(motif_name) = NULL
  colnames(motif_pk_match) = motif_name      
}

saveRDS(motif_ix, file = 'TRN/Motifs/BALL_motif_match_pvE-03.rds')
saveRDS(motif_ix_str, file = 'TRN/Motifs/BALL_motif_match_pv5E-05.rds')


## check motif overlapping with any peak -- T-ALL ####

seurat.atac.TALL <- readRDS(file = 'PanLeukemia/scATAC/Motif/AllTALL_scATAC_downsampleBySeuratProjection_Signca_RecalledPeaks_BlacklistRemoved_Oct12_2022_withChromVAR.rds')

DefaultAssay(seurat.atac.TALL) <- 'ATAC'

peaks = rownames(seurat.atac.TALL)
peaks = sapply(peaks, function(x) unlist(strsplit(x, ','))[1])
names(peaks) = NULL

pks = tidyr::separate(data = data.table('peak_name' = peaks),
                      col = 'peak_name', into = c('chr', 'start', 'end'),
                      remove = F)
pks$start = as.integer(pks$start)
pks$end = as.integer(pks$end)
setkey(pks, chr, start)
# Make a set of peaks
peaks <- GenomicRanges::GRanges(seqnames = pks$chr,
                                ranges = IRanges::IRanges(start = pks$start,
                                                          end = pks$end))

motif_ix <- matchMotifs(human_pwms_v2, peaks, 
                        genome = BSgenome.Hsapiens.UCSC.hg38,
                        p.cutoff = 1e-3)
motif_ix_str <- matchMotifs(human_pwms_v2, peaks, 
                            genome = BSgenome.Hsapiens.UCSC.hg38,
                            p.cutoff = 5e-5)

if(F){
  motif_name = motif_ix@colData$name
  pk_inf = motif_ix@rowRanges
  motif_pk_match <- motif_ix@assays@data$motifMatches
  rownames(motif_pk_match) = paste0(pk_inf@seqnames, '-', pk_inf@ranges)
  names(motif_name) = NULL
  colnames(motif_pk_match) = motif_name      
}

saveRDS(motif_ix, file = 'TRN/Motifs/AllTALL_motif_match_pvE-03.rds')
saveRDS(motif_ix_str, file = 'TRN/Motifs/AllTALL_motif_match_pv5E-05.rds')

