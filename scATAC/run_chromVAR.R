source('dsAnalysis_utilities.R')
library(Seurat)
library(Signac)
library(motifmatchr)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)
`%notin%` = Negate(`%in%`)

##reture chromvar result given a seurat obj

seuratFile <- 'PanLeukemia/scATAC/SeuratObj/scATAC_Signca_RecalledPeaks_BlacklistRemoved.rds'
seurat.atac <- readRDS(seuratFile)

mtx <- seurat.atac@assays$ATAC@counts

chr.exclude.1 <- rownames(mtx)[grep("chrM", rownames(mtx))]
chr.exclude.2 <- rownames(mtx)[grep("chrUn", rownames(mtx))]
chr.exclude.3 <- rownames(mtx)[grep("random", rownames(mtx))]

chr.exclude <- c(chr.exclude.1, chr.exclude.2, chr.exclude.3)

peak.include <- setdiff(rownames(mtx), chr.exclude)
mtx <- mtx[peak.include,]

## further filter peaks
rs <- Matrix::rowSums(mtx > 0)
mtx <- mtx[rs > 0.005, ]

chromVAR.obj = run_chromVAR(mtx, genomeName = 'BSgenome.Hsapiens.UCSC.hg38')

chromvarPath <- 'PanLeukemia/scATAC/chromVAR_Obj/ChromVAR_Related2_scATAC_Signca_RecalledPeaks_BlacklistRemoved.rds'
saveRDS(chromVAR.obj, file = chromvarPath)

