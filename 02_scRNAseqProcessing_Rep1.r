#Load libraries
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(reshape2)
library(pheatmap)
library(proxyC)
library(Matrix)
library(igraph)
library(DoubletFinder)
library(harmony)
library(reshape2)
library(dplyr)
library(SeuratDisk)

#Convert h5ad to seurat
Convert("./197397/counts_filtered/adata.h5ad", dest = "h5seurat", overwrite = TRUE) 
seuratObject <- LoadH5Seurat("./197397/counts_filtered/adata.h5seurat")


#Filter cells
##Add mitochondrial content
seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, pattern = "^AMEXmito",assay = "spliced")

##Apply cells
seuratObject <- subset(x = seuratObject, subset = nCount_spliced > 5000 & nFeature_spliced >1000 & percent.mt < 5)

#Log normalize the data
seuratObject <- NormalizeData(seuratObject, normalization.method = "LogNormalize", scale.factor = 10000)

#Find the 5000 most variable genes
seuratObject <- FindVariableFeatures(seuratObject, selection.method = "vst", nfeatures = 5000)

#Scale the data so that the mean of each gene is 0
all.genes <- rownames(seuratObject)
seuratObject <- ScaleData(seuratObject, features = all.genes)

#Find Principal components of the dataset
seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject),verbose = F)

#Find neighbors for UMAP and clusters for all the cells
seuratObject <- FindNeighbors(seuratObject, dims = 1:10,verbose = FALSE)
seuratObject <- FindClusters(seuratObject, resolution = 0.5,verbose = FALSE)

#Perform dimensionality reduction using UMAP
seuratObject <- RunUMAP(seuratObject, dims = 1:10,verbose = F)

#Plot UMAP
DimPlot(seuratObject, reduction = "umap")

#Regress cell cycle effects
##Score cell cycle
s.genes <- rownames(Annotation[Annotation[,1] %in% cc.genes$s.genes,])
g2m.genes <- rownames(Annotation[Annotation[,1] %in% cc.genes$g2m.genes,])
seuratObject <- CellCycleScoring(seuratObject, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

##Regress cell cycle
seuratObject <- ScaleData(seuratObject, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seuratObject))

##Re-do UMAP
seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject))
seuratObject <- FindNeighbors(seuratObject, dims = 1:10)
seuratObject <- FindClusters(seuratObject, resolution = 0.5)
seuratObject <- RunUMAP(seuratObject, dims = 1:10)

#Remove doublets
##Perform sweep
Sweep_seurat <- paramSweep_v3(seuratObject, PCs = 1:10, sct = FALSE)
sweep.stats_seurat <- summarizeSweep(Sweep_seurat, GT = FALSE)
pk_seurat <- find.pK(sweep.stats_seurat)

##Find doublets
seurat_annotations    <- seuratObject$seurat_clusters
homotypic.prop_seurat <- modelHomotypic(seurat_annotations)
nExp_poi              <- round(0.075*nrow(seuratObject@meta.data))
nExp_poi.adj          <- round(nExp_poi*(1-homotypic.prop_seurat))

seuratObject_Doublets<- doubletFinder_v3(seuratObject, PCs = 1:10, pN = 0.25, pK = 0.04, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

##Remove doublets
seuratObject <- seuratObject[,seuratObject_Doublets$DF.classifications_0.25_0.04_585 == "Singlet"]

#Save Object
saveRDS(seuratObject,file="./SeuratRep1.RDS")