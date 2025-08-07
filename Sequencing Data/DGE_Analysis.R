# load required libraries
library(Seurat)
library(tidyverse)
library(DESeq2)


## SMA vs control

### Angelica samples
sp = c('OC63','OC80','OC287','OC111','OC119','OC122')
org = 'd30NMO'
ct = 'Neurons'

# Filter Seurat object
dat <- subset(x = mydata, subset = sample %in% sp & cell.types2 %in% ct)
dat$sample2 <- ifelse(dat$sample %in% c('OC63','OC80'), 'OC63.80', dat$sample)
table(as.character(dat$cell.types2))
table(dat$sample,dat$subject);table(dat$sample,dat$group);table(dat$sample2,dat$sample)

# filter out genes with low expression in all samples
min.expr.rate = 0.1
counts <- SeuratObject::LayerData(dat,assay='RNA', layer='counts')
low.genes <- lapply(unique(dat$sample),function(i){
  my.counts <- counts[,grep(x= colnames(counts), pattern = i, value=T)]
  my.genes <- rownames(my.counts[rowMeans(my.counts > 0) < min.expr.rate,])
  my.genes
})
selected.genes <- setdiff(rownames(dat), purrr::reduce(low.genes, intersect))

# pseudobulk the counts
pb <- AggregateExpression(dat, assays = "RNA",features = selected.genes, return.seurat = F, group.by = c("sample2", "group"))
pb <- as.matrix(pb$RNA);dim(pb)

md <- lapply(str_split(string = colnames(pb), pattern = '_'), function(x) data.frame(sample2= x[1], group = x[2])) %>% list_rbind
md

# define reference level of the group variable
md$group <- factor(md$group, levels = c('ctrl','sma'))
#md$sample <- factor(md$sample, levels = sp)
md$sample2 <- factor(md$sample2, levels = c('OC63.80','OC287','OC111','OC119','OC122'))
table(md$group);table(md$sample2)

# run DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = pb,
                              colData = md,
                              design = ~ group)

# dds <- DESeqDataSetFromMatrix(countData = pb,
#                               colData = md,
#                               design = ~ sample)
# This does not work because dispersion can not be estimated
dds
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds,name = "group_sma_vs_ctrl") %>% as.data.frame()
res2 <- lfcShrink(dds, coef="group_sma_vs_ctrl", type="apeglm") %>% as.data.frame()

# save results 
saveRDS(res2, paste(org,'Neurons','Angelica','deseq2_sma_vs_ctrl.rds', sep = '_'))

### Ines samples
dat <- subset(x = mydata, subset = cell.types2 %in% ct & sample %in% c('O48','O51'))
#dat <- subset(x = mydata, subset = cell.types2 %in% 'Neuron' & sample %in% c('O48','O51'))
table(dat$sample);table(as.character(dat$cell.types2))
# filter out genes with low expression in all samples
min.expr.rate = 0.1
counts <- SeuratObject::LayerData(dat,assay='RNA', layer='counts')
low.genes <- lapply(unique(dat$sample),function(i){
  my.counts <- counts[,grep(x= colnames(counts), pattern = i, value=T)]
  my.genes <- rownames(my.counts[rowMeans(my.counts > 0) < min.expr.rate,])
  my.genes
})
selected.genes <- setdiff(rownames(dat), purrr::reduce(low.genes, intersect))

Idents(dat) <- 'sample'
res3 <- FindMarkers(dat, ident.1 = "O51", ident.2 = "O48", features = selected.genes, logfc.threshold = 0, min.pct = 0)

# save results 
saveRDS(res3, paste(org,'Neurons','Ines','wilcoxon_sma_vs_ctrl.rds', sep = '_'))
