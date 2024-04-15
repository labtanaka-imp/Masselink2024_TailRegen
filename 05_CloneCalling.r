library(Seurat)
library(ggplot2)
library(reshape2)
library(proxyC)
library(Matrix)
library(harmony)
library(reshape2)
library(dplyr)

# Load both Matrixes and combine them
seuratRep1<- readRDS("./SeuratRep1.RDS")
seuratRep2  <- readRDS("./SeuratRep2.RDS")
seuratCombo <- merge(seuratRep1,y=seuratRep2,add.cell.ids = c("Rep1","Rep2"))

#Perform normalization and calculate UMAP
##Log normalize the data
seuratCombo <- NormalizeData(seuratCombo, normalization.method = "LogNormalize", scale.factor = 10000)

##Find the 5000 most variable genes
seuratCombo <- FindVariableFeatures(seuratCombo, selection.method = "vst", nfeatures = 5000)

##Scale the data so that the mean of each gene is 0
all.genes <- rownames(seuratCombo)
seuratCombo <- ScaleData(seuratCombo, features = all.genes)

##Find Principal components of the dataset
seuratCombo <- RunPCA(seuratCombo, features = VariableFeatures(object = seuratCombo),verbose = F)

#F#ind neighbors for UMAP and clusters for all the cells
seuratCombo <- FindNeighbors(seuratCombo, dims = 1:10,verbose = FALSE)
seuratCombo <- FindClusters(seuratCombo, resolution = 0.5,verbose = FALSE)

##Perform dimensionality reduction using UMAP
seuratCombo <- RunUMAP(seuratCombo, dims = 1:10,verbose = F)

#Harmonize the data
seuratCombo <- RunHarmony(seuratCombo, "dataset")
seuratCombo <- RunUMAP(seuratCombo, reduction = "harmony",dims = 1:20)



###########Barcodes
# Load Barcodes into memory
##Set Barcodes
BarcodesV <- c('ACCGGT','CATCAC','CGTACA')
Versions  <- c('V1','V2','V3')

##All were filtered by selecting those with at least two counts per cell
Barcodes_Rep1_1 <- read.table('./Rep1.Uniq.FilteredCounts.txt')
Barcodes_Rep2_1 <- read.table('./Rep2_1.Uniq.FilteredCounts.txt')
Barcodes_Rep2_2 <- read.table('./Rep2_2.Uniq.FilteredCounts.txt')

##Create Binary matrixes
MatrixRep1 <- list()
MatrixRep2_1 <- list()
MatrixRep2_2 <- list()

for (i in 1:3){
    MatrixRep1[[Versions[i]]] <- matrix(0,length(unique(Barcodes_Rep1_1[Barcodes_Rep1_1[,3] == BarcodesV[i],2])),length(unique(Barcodes_Rep1_1[Barcodes_Rep1_1[,3] == BarcodesV[i],4])))
    rownames(MatrixRep1[[Versions[i]]]) <- unique(Barcodes_Rep1_1[Barcodes_Rep1_1[,3] == BarcodesV[i],2])
    colnames(MatrixRep1[[Versions[i]]]) <- unique(Barcodes_Rep1_1[Barcodes_Rep1_1[,3] == BarcodesV[i],4])
}

for (i in 1:3){
    MatrixRep2_1[[Versions[i]]] <- matrix(0,length(unique(Barcodes_Rep2_1[Barcodes_Rep2_1[,3] == BarcodesV[i],2])),length(unique(Barcodes_Rep2_1[Barcodes_Rep2_1[,3] == BarcodesV[i],4])))
    rownames(MatrixRep2_1[[Versions[i]]]) <- unique(Barcodes_Rep2_1[Barcodes_Rep2_1[,3] == BarcodesV[i],2])
    colnames(MatrixRep2_1[[Versions[i]]]) <- unique(Barcodes_Rep2_1[Barcodes_Rep2_1[,3] == BarcodesV[i],4])
}

for (i in 1:3){
    MatrixRep2_2[[Versions[i]]] <- matrix(0,length(unique(Barcodes_Rep2_2[Barcodes_Rep2_2[,3] == BarcodesV[i],2])),length(unique(Barcodes_Rep2_2[Barcodes_Rep2_2[,3] == BarcodesV[i],4])))
    rownames(MatrixRep2_2[[Versions[i]]]) <- unique(Barcodes_Rep2_2[Barcodes_Rep2_2[,3] == BarcodesV[i],2])
    colnames(MatrixRep2_2[[Versions[i]]]) <- unique(Barcodes_Rep2_2[Barcodes_Rep2_2[,3] == BarcodesV[i],4])
}


##Fill where there is barcodes
for (j in 1:3){
    for (i in unique(Barcodes_Rep1_1[Barcodes_Rep1_1[,3] == BarcodesV[j],4])){
         MatrixRep1[[Versions[j]]][Barcodes_Rep1_1[(Barcodes_Rep1_1[,4]==i) & (Barcodes_Rep1_1[,3]==BarcodesV[j]),2],i] <- 1
    }
}

for (j in 1:3){
    for (i in unique(Barcodes_Rep2_1[Barcodes_Rep2_1[,3] == BarcodesV[j],4])){
         MatrixRep2_1[[Versions[j]]][Barcodes_Rep2_1[(Barcodes_Rep2_1[,4]==i) & (Barcodes_Rep2_1[,3]==BarcodesV[j]),2],i] <- 1
    }
}

for (j in 1:3){
    for (i in unique(Barcodes_Rep2_2[Barcodes_Rep2_2[,3] == BarcodesV[j],4])){
         MatrixRep2_2[[Versions[j]]][Barcodes_Rep2_2[(Barcodes_Rep2_2[,4]==i) & (Barcodes_Rep2_2[,3]==BarcodesV[j]),2],i] <- 1
    }
}

##Create list where results will be stored
ProcessedMatrix_Rep1 <- list()
ProcessedMatrix_Rep2 <- list()
Clones_Rep1 <- list()
Clones_Rep2 <- list()

##Rep1
for (i in Versions) {
    CurrVer <- Versions[i]
    print(dim(MatrixRep1[[CurrVer]]))

    #Remove cells that arenot in the dataset
    MatrixRep1[[CurrVer]] <- MatrixRep1[[CurrVer]][,colnames(MatrixRep1[[CurrVer]]) %in% gsub(pattern = "Rep1_",replacement = "",Cells(seuratCombo_Reduced)[grep(pattern = "Rep1_",Cells(seuratCombo_Reduced))])]
    print(dim(MatrixRep1[[CurrVer]]))

    #Remove barcodes in only one cells
    MatrixRep1[[CurrVer]] <- MatrixRep1[[CurrVer]][rowSums(MatrixRep1[[CurrVer]]) > 1,]
    print(dim(MatrixRep1[[CurrVer]]))

    #Remove 90th quantile
    print(quantile(rowSums(MatrixRep1[[CurrVer]]),0.90))
    MatrixRep1[[CurrVer]] <- MatrixRep1[[CurrVer]][rowSums(MatrixRep1[[CurrVer]]) < quantile(rowSums(MatrixRep1[[CurrVer]]),0.90),]
    print(dim(MatrixRep1[[CurrVer]]))

    #Remove barcodes in only one cells
    MatrixRep1[[CurrVer]] <- MatrixRep1[[CurrVer]][rowSums(MatrixRep1[[CurrVer]]) > 1,]
    print(dim(MatrixRep1[[CurrVer]]))

    #Remove Empty Cells and with more than 20 barcodes
    MatrixRep1[[CurrVer]] <- MatrixRep1[[CurrVer]][,(colSums(MatrixRep1[[CurrVer]]) > 1) & (colSums(MatrixRep1[[CurrVer]]) <= 20)]
    print(dim(MatrixRep1[[CurrVer]]))

    #Remove barcodes in only one cells
    MatrixRep1[[CurrVer]] <- MatrixRep1[[CurrVer]][rowSums(MatrixRep1[[CurrVer]]) > 4,]
    MatrixRep1[[CurrVer]] <- MatrixRep1[[CurrVer]][,(colSums(MatrixRep1[[CurrVer]]) > 1) ]
    print(dim(MatrixRep1[[CurrVer]]))
    ProcessedMatrix_Rep1[[CurrVer]] <- MatrixRep1[[CurrVer]]

    #Calculate Jaccard Index
    Jac <- simil(t(MatrixRep1[[CurrVer]]),method = "jaccard")

    jac.summ <- Matrix::summary(Jac)
    jac.lower.i <- jac.summ$j
    jac.summ$j <- jac.summ$i
    jac.summ$i <- jac.lower.i
    lower.tri.summ <- subset(jac.summ, i>j) # Exclude diagnol


    Jac_Matrix <- sparseMatrix(i = lower.tri.summ$i,
                            j = lower.tri.summ$j,
                            x = lower.tri.summ$x,
                            dims = dim(Jac),
                            dimnames = dimnames(Jac))

    Jac.df     <- as.data.frame(Matrix::summary(Jac_Matrix))
    Jac.df.sub <- Jac.df[which(Jac.df$x > 0.45), ]
    check.corelation <- Jac.df.sub[,c(1,2)]
    colnames(check.corelation) <- c("row", "col")
    check.corelation <- as.matrix(check.corelation)

    graph.cor <- graph.data.frame(check.corelation, directed = FALSE)
    groups.cor <- split(unique(as.vector(check.corelation)), clusters(graph.cor)$membership)
    Clones_Rep1[CurrVer] <- lapply(groups.cor,
                            function(list.cor){
                            rownames(Jac_Matrix)[list.cor]})
}

##Rep2
for (i in Versions) {
    CurrVer <- Versions[i]
    print(dim(MatrixRep2_2[[CurrVer]]))

    #Remove cells that arenot in the dataset
    MatrixRep2_2[[CurrVer]] <- MatrixRep2_2[[CurrVer]][,colnames(MatrixRep2_2[[CurrVer]]) %in% gsub(pattern = "Rep2_",replacement = "",Cells(seuratCombo_Reduced)[grep(pattern = "Rep2_",Cells(seuratCombo_Reduced))])]
    print(dim(MatrixRep2_2[[CurrVer]]))

    #Remove barcodes in only one cells
    MatrixRep2_2[[CurrVer]] <- MatrixRep2_2[[CurrVer]][rowSums(MatrixRep2_2[[CurrVer]]) > 1,]
    print(dim(MatrixRep2_2[[CurrVer]]))

    #Remove 90th quantile
    print(quantile(rowSums(MatrixRep2_2[[CurrVer]]),0.90))
    MatrixRep2_2[[CurrVer]] <- MatrixRep2_2[[CurrVer]][rowSums(MatrixRep2_2[[CurrVer]]) < quantile(rowSums(MatrixRep2_2[[CurrVer]]),0.90),]
    print(dim(MatrixRep2_2[[CurrVer]]))

    #Remove barcodes in only one cells
    MatrixRep2_2[[CurrVer]] <- MatrixRep2_2[[CurrVer]][rowSums(MatrixRep2_2[[CurrVer]]) > 1,]
    print(dim(MatrixRep2_2[[CurrVer]]))

    #Remove Empty Cells and with more than 20 barcodes
    MatrixRep2_2[[CurrVer]] <- MatrixRep2_2[[CurrVer]][,(colSums(MatrixRep2_2[[CurrVer]]) > 1) & (colSums(MatrixRep2_2[[CurrVer]]) <= 20)]
    print(dim(MatrixRep2_2[[CurrVer]]))

    #Remove barcodes in only one cells
    MatrixRep2_2[[CurrVer]] <- MatrixRep2_2[[CurrVer]][rowSums(MatrixRep2_2[[CurrVer]]) > 4,]
    MatrixRep2_2[[CurrVer]] <- MatrixRep2_2[[CurrVer]][,(colSums(MatrixRep2_2[[CurrVer]]) > 1) ]
    print(dim(MatrixRep2_2[[CurrVer]]))
    ProcessedMatrix_Rep2_2[[CurrVer]] <- MatrixRep2_2[[CurrVer]]

    #Calculate Jaccard Index
    Jac <- simil(t(MatrixRep2_2[[CurrVer]]),method = "jaccard")

    jac.summ <- Matrix::summary(Jac)
    jac.lower.i <- jac.summ$j
    jac.summ$j <- jac.summ$i
    jac.summ$i <- jac.lower.i
    lower.tri.summ <- subset(jac.summ, i>j) # Exclude diagnol


    Jac_Matrix <- sparseMatrix(i = lower.tri.summ$i,
                            j = lower.tri.summ$j,
                            x = lower.tri.summ$x,
                            dims = dim(Jac),
                            dimnames = dimnames(Jac))

    Jac.df     <- as.data.frame(Matrix::summary(Jac_Matrix))
    Jac.df.sub <- Jac.df[which(Jac.df$x > 0.45), ]
    check.corelation <- Jac.df.sub[,c(1,2)]
    colnames(check.corelation) <- c("row", "col")
    check.corelation <- as.matrix(check.corelation)

    graph.cor <- graph.data.frame(check.corelation, directed = FALSE)
    groups.cor <- split(unique(as.vector(check.corelation)), clusters(graph.cor)$membership)
    Clones_Rep2_2[CurrVer] <- lapply(groups.cor,
                            function(list.cor){
                            rownames(Jac_Matrix)[list.cor]})
}