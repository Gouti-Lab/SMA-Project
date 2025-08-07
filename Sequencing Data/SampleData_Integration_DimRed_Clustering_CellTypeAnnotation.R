# load required libraries
library(tidyverse)
library(Seurat)

# set globals maxSize
options(future.globals.maxSize = 16 * 1024^3)

# read the day 30 NMO snRNA-seq filtered dataset
fil.dataset <- readRDS('d30NMO_snRNAseq_filtered_dataset.rds')

# Please note that our full dataset includes more NMO samples from other lines and timepoints, 
# but we have only included here day 30 samples that were used in the current manuscript.
# This is why running PCA followed by UMAP with this subset of samples may result in a different UMAP embedding shape that the one presented in the manuscript.
# However, this will not change the results of cell type annotations or any other downstream analysis carried out. 

# run Seurat pipeline
# normalizing the data
fil.dataset <- NormalizeData(fil.dataset, normalization.method = "LogNormalize", scale.factor = 10000)

# identify highly variable features
fil.dataset <- FindVariableFeatures(fil.dataset, selection.method = "vst", nfeatures = 5000)

# scaling the data
fil.dataset <- ScaleData(fil.dataset)

# Run linear dimensional reduction
fil.dataset <- RunPCA(fil.dataset)

# cluster cells
res = seq(0.4,1.2,by = 0.1) # resolutions to test 
fil.dataset <- FindNeighbors(fil.dataset, dims = 1:50, reduction = "pca")
fil.dataset <- FindClusters(fil.dataset, resolution = res, cluster.name = paste("unintegrated_clusters",res,sep = '_'))

# Run non-linear dimensional reduction
fil.dataset <- RunUMAP(fil.dataset, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")

# Perform Harmony integration
fil.dataset <- IntegrateLayers(
  object = fil.dataset, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE
)

# cluster cells after integration
fil.dataset <- FindNeighbors(fil.dataset, reduction = "harmony", dims = 1:50)
fil.dataset <- FindClusters(fil.dataset, resolution = res, cluster.name = paste("harmony_clusters",res,sep = '_'))

# Run non-linear dimensional reduction after integration
fil.dataset <- RunUMAP(fil.dataset, reduction = "harmony", dims = 1:50, reduction.name = "umap.harmony")

# Perform CCA integration
fil.dataset <- IntegrateLayers(
  object = fil.dataset, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

# cluster cells after integration
fil.dataset <- FindNeighbors(fil.dataset, reduction = "integrated.cca", dims = 1:50)
fil.dataset <- FindClusters(fil.dataset, resolution = res, cluster.name = paste("cca_clusters",res,sep = '_'))

# Run non-linear dimensional reduction after integration
fil.dataset <- RunUMAP(fil.dataset, reduction = "integrated.cca", dims = 1:50, reduction.name = "umap.cca")

# Perform RPCA integration
fil.dataset <- IntegrateLayers(
  object = fil.dataset, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

# cluster cells after integration
fil.dataset <- FindNeighbors(fil.dataset, reduction = "integrated.rpca", dims = 1:50)
fil.dataset <- FindClusters(fil.dataset, resolution = res, cluster.name = paste("rpca_clusters",res,sep = '_'))

# Run non-linear dimensional reduction after integration
fil.dataset <- RunUMAP(fil.dataset, reduction = "integrated.rpca", dims = 1:50, reduction.name = "umap.rpca")

# Join the layers of the seurat object following integration
fil.dataset <- JoinLayers(fil.dataset)

# save the integrated object 
saveRDS(fil.dataset, 'd30NMO_snRNAseq_integrated_dataset.rds')

# After testing different integration algorithms and clusters obtained at different resolutions, Harmony clusters at a resolution of 0.5 were selected for cell type annotations

# To enhance separation of clusters in the UMAP embedding, UMAP was re-run with different parameters (optional step)
fil.dataset <- RunUMAP(fil.dataset, reduction = "harmony", dims = 1:50, reduction.name = "umap.harmony2", min.dist = 0.2, n.neighbors = 50)

# set identities of the seurat object before testing for cluster biomarkers
Idents(fil.dataset) <- 'harmony_clusters_0.5'

# Cluster biomarker identification
## find markers for every cluster compared to all remaining cells, report only the positive ones with log2FC > 0.1 and p-value < 0.01
cluster.markers <- FindAllMarkers(fil.dataset, only.pos = TRUE, logfc.threshold = 0.1, return.thresh = 0.01)

# Cell type annotation of the clusters based on the identified biomarkers
fil.dataset$cell.types <- case_match(fil.dataset$`harmony_clusters_0.5`, '0' ~ 'Neurons 1', '1' ~ 'Low RNA 1', '2' ~ 'Fibroblasts 1', '3' ~ 'Glia', '4' ~ 'Neurons 2', '5' ~ 'Ependymal cells', '6' ~ 'Immature neurons 1', 
                                '7' ~ 'Satellite cells', '8' ~ 'Neurons 3', '9' ~ 'Neural progenitors', '10' ~ 'Neurons 4', '11' ~ 'Skeletal muscle fibers', '12' ~ 'Unidentifiable 1', '13' ~ 'Low RNA 2', '14' ~ 'Myocytes',
                                '15' ~ 'Fibroblasts 2', '16' ~ 'Proliferating fibroblasts', '17' ~ 'Myogenic progenitors', '18' ~ 'Unidentifiable 2', '19' ~ 'Immature neurons 2', '20' ~ 'Fibroblasts 3', '21' ~ 'Low RNA 3',
                                '22' ~ 'Epithelial cells', '23' ~ 'Neurons 5', '24' ~ 'Podocytes', '25' ~ 'Schwann cells','26' ~ 'Endothelial cells')

fil.dataset$cell.types <- factor(fil.dataset$cell.types, levels = c('Neural progenitors', paste('Immature neurons',1:2,sep = ' '),
                                                          paste('Neurons',1:5, sep = ' '), 'Glia','Schwann cells','Ependymal cells',
                                                          'Myogenic progenitors','Satellite cells','Myocytes','Skeletal muscle fibers',
                                                          'Proliferating fibroblasts',paste('Fibroblasts',1:3,sep = ' '),
                                                          'Epithelial cells', 'Endothelial cells','Podocytes',
                                                          paste('Low RNA',1:3,sep = ' '),
                                                          paste('Unidentifiable',1:2,sep = ' '))) 

# Aggregate similar clusters 
fil.dataset$cell.types2 <- case_match(fil.dataset$`harmony_clusters_0.5`, '0' ~ 'Neurons', '1' ~ 'Low RNA', '2' ~ 'Fibroblasts', '3' ~ 'Glia', '4' ~ 'Neurons', '5' ~ 'Ependymal cells', '6' ~ 'Immature neurons', 
                                 '7' ~ 'Satellite cells', '8' ~ 'Neurons', '9' ~ 'Neural progenitors', '10' ~ 'Neurons', '11' ~ 'Skeletal muscle fibers', '12' ~ 'Unidentifiable', '13' ~ 'Low RNA', '14' ~ 'Myocytes',
                                 '15' ~ 'Fibroblasts', '16' ~ 'Proliferating fibroblasts', '17' ~ 'Myogenic progenitors', '18' ~ 'Unidentifiable', '19' ~ 'Immature neurons', '20' ~ 'Fibroblasts', '21' ~ 'Low RNA',
                                 '22' ~ 'Epithelial cells', '23' ~ 'Neurons', '24' ~ 'Podocytes', '25' ~ 'Schwann cells','26' ~ 'Endothelial cells')

fil.dataset$cell.types2 <- factor(fil.dataset$cell.types2, levels = c('Neural progenitors',
                                                            'Immature neurons',
                                                            'Neurons', 'Glia','Schwann cells','Ependymal cells',
                                                            'Myogenic progenitors','Satellite cells','Myocytes','Skeletal muscle fibers',
                                                            'Proliferating fibroblasts','Fibroblasts',
                                                            'Epithelial cells', 'Endothelial cells','Podocytes',
                                                            'Low RNA',
                                                            'Unidentifiable')) 


# save the integrated and annotated object 
saveRDS(fil.dataset, 'd30NMO_snRNAseq_integrated_annotated_dataset.rds')

