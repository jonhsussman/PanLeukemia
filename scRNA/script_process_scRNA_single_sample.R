source('scDataAnalysis_Utilities.R')


## path of cellrangerange result for  each sample ##
args = commandArgs(T)
dir0 = "/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/CellRangerResults/June26_2020/"
sampleName = args[1]
## process one by one
path2mtx = paste0(dir0, sampleName, '_scRNA/outs/filtered_feature_bc_matrix')
seurat.obj = doSeurat_rmDoublets_dir(filtered_mtx_dir = path2mtx, 
                                     0.05, pmito.upper = 0.15, 
                                     min_nCount = 1500,
                                     max_nCount = 40000,
                                     bcPrefix = sampleName)


saveRDS(seurat.obj, file = paste0('Seurat_Objects/scRNA/HealthyDonor/seurat_', sampleName, '_doubletRemoved.rds'))


