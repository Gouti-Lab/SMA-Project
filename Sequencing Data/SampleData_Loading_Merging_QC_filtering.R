# load necessary packages
library(tidyverse)
library(Seurat)

# define sample IDs
samples <- c('O48','O51','O54','O57','O72','O75','O78','O81','O84','O87','OC63','OC80','OC287','OC111','OC119','OC122')

# Please note that our full dataset includes more NMO samples from other lines and timepoints, 
# but we have only included here day 30 samples that were used in the current manuscript.
# This is why running PCA followed by UMAP with this subset of samples may result in a different UMAP embedding shape that the one presented in the manuscript.
# However, this will not change the results of cell type annotations or any other downstream analysis carried out. 

# save the filtered feature-barcode matrices of all samples as a list of seurat objects
dir <- '/path/to/CellRanger_Output/'
obj <- list()
for (i in samples) {
  # read in data from the Cellranger output folder
  data <- Read10X(paste0(dir,i,'/filtered_feature_bc_matrix/'))
  # create seurat object from data, filter out nuclei with less than 100 features (genes)
  obj[[i]] <- CreateSeuratObject(counts = data, min.features = 100)
  # add sample annotation to the object
  obj[[i]]$sample <- i
} 

# merge the list into a single seurat object
merged.obj <- merge(obj[[1]], obj[2:length(obj)], add.cell.ids= names(obj)) # merge with addition of sample IDs to cell IDs to make them unique

# add the organoid line annotations to the merged object metadata
merged.obj$subject <- case_match(merged.obj$sample,
                                 c('OC63','OC80','O48') ~ 'ctrl1',
                                 'OC287' ~ 'ctrl2',
                                 c('OC111','O51','O54','O57','O72','O75','O78','O81','O84','O87') ~ 'sma1',
                                 'OC119' ~ 'sma2',
                                 'OC122' ~ 'sma3')
# add timepoint annotation to the merged object metadata
merged.obj$timepoint <- 'd30'

# add treatment annotations to the merged object metadata
merged.obj$treatment <- case_match(merged.obj$sample,
                                   c('O57','O75','O84') ~ 'risdiplam_250nM',
                                   c('O54','O78','O87') ~ 'branaplam_40nM',
                                   .default = 'none'
)

# add researcher annotations
merged.obj$researcher <- case_match(merged.obj$sample,
                                    c('O48','O51','O54','O57','O72','O75','O78','O81','O84','O87') ~ 'Ines',
                                    .default = 'Angelica'
)

# calculate percentage of mitochondrial reads for each cell/nucleus
merged.obj[['percent.mt']] <- PercentageFeatureSet(merged.obj, pattern = 'MT-')

# exclude low-quality nuclei
## Nuclei with a percentage of mitochondrial reads greater than 2% and a feature count less than 100 were considered of low quality and excluded
## Nuclei with a feature count greater than 6000 were considered doublets or aggregates and excluded
fil.dataset <- subset(merged.obj, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt <= 2)

# check percentage of nuclei excluded after QC filtering
100*( ncol(merged.obj) - ncol(fil.dataset) ) / ncol(merged.obj) # almost 6% of nuclei were excluded at this stage

# save filtered seurat object after QC
saveRDS(fil.dataset, 'd30NMO_snRNAseq_filtered_dataset.rds')

