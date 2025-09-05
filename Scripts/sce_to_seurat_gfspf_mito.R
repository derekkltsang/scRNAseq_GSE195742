#-----------------------------Install required Pkgs-----------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Seurat")
BiocManager::install("dplyr")
BiocManager::install("Matrix")
BiocManager::install("scater")
BiocManager::install("scran")
BiocManager::install("ggraph")
BiocManager::install("ggpubr")

install.packages("remotes")
remotes::install_github(repo ='mojaveazure/loomR', ref = 'develop')

library(scran)
library(dplyr)
library(scater)
library(Seurat)
library(Matrix)
library(ggplot2)
library(loomR)
library(magrittr)
library(ggpubr)
library(org.Mm.eg.db)
library(cowplot)

#----------------------------- Generate Seurat Object -----------------------------------

#  Make Seurat object
sce.data <- Read10X("/home/derek/Documents/R_Projects/CR_SCR_raw_matrix/S02_FM")
seurat.sce  <- new("seurat", raw.data = sce.data)
seurat.sce <- CreateSeuratObject(counts = sce.data, project = "S02")
#seurat.sce <- CreateSeuratObject(counts = sce.data, project = "SPF_WT", min.cells = 3, min.features = 200)

# Make SCE Object
sce <- as.SingleCellExperiment(seurat.sce)
#m_srt <- readRDS(file = "S02_QC.rds")

# Low gene expression removal (genes expressed in only 5 cells)
numcells <- nexprs(sce, byrow=TRUE)
keep <- numcells >= 5
sum(keep)

sce <- sce[keep]

hist(log10(numcells), breaks=100, main="", col="grey80", 
     xlab=expression("average count"))
abline(v=log10(40), col="blue", lwd=2, lty=2)

# ID cells with high mitochondria gene expression
is.mito <- grepl("^mt-", rownames(sce))
#?perCellQCMetrics
#sce <- perCellQCMetrics(sce)
sce <- addPerCellQC(sce, subsets=list(Mito=grep("mt-", rownames(sce)))) 

sce$subsets_Mito_percent

CtThres = 0
FeatThres = 0
MtThres = 100

rem.total <- sce$nCount_RNA < CtThres
rem.n <- sce$nFeature_RNA < FeatThres
rem.mt <- sce$subsets_Mito_percent > MtThres
sum(rem.total)
sum(rem.n)
sum(rem.mt)

#Make QC Plots
png("S02raw_mito_QCplot.png", width =2000, height =500, units = "px", pointsize = 17)
par(mfrow=c(1,4))
#Count
hist(sce$nCount_RNA, xlab="Library size", main="",
     breaks=50, col="grey80", ylab="Number of cells")
abline(v=CtThres, col="blue", lwd=2, lty=2)

hist(sce$nCount_RNA, xlab="Library size", main="",
     breaks=800, col="grey80", xlim=c(0,5000), ylab="Number of cells")
abline(v=CtThres, col="blue", lwd=2, lty=2)
#Features
hist(sce$nFeature_RNA, xlab="Number of expressed genes", main="",
     breaks=100, col="grey80", xlim=c(0,5000), ylab="Number of cells")
abline(v=FeatThres, col="blue", lwd=2, lty=2)
#Mito
hist(sce$subsets_Mito_percent, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=40, main="", col="grey80")
abline(v=MtThres, col="blue", lwd=2, lty=2)
dev.off()

p4 <- plotColData(sce, x="nCount_RNA", y="nFeature_RNA") +
  geom_vline(xintercept=CtThres, color="red", size=1, linetype="dashed") +
  geom_hline(yintercept=FeatThres, color="red", size=1, linetype="dashed")

p5 <- plotColData(sce, x="nCount_RNA", y="subsets_Mito_percent") + 
  geom_vline(xintercept=CtThres, color="red", size=1, linetype="dashed") +
  geom_hline(yintercept=MtThres, color="red", size=1, linetype="dashed")

p6 <- plotColData(sce, x="nFeature_RNA", y="subsets_Mito_percent") +
  geom_vline(xintercept=FeatThres, color="red", size=1, linetype="dashed") +
  geom_hline(yintercept=MtThres, color="red", size=1, linetype="dashed")

png("S02raw_mito_QCplots2.png", width =1500, height =500, units = "px", pointsize = 17)
plot_grid(p4,p5,p6,
          nrow=1, ncol=3)
dev.off()

sce <- sce[,!(rem.total | rem.n | rem.mt)]

filter.sce <- data.frame(ByLibSize=sum(rem.total), ByFeature=sum(rem.n), ByMito=sum(rem.mt), Remaining=ncol(sce),
                         CountThres=as.numeric(CtThres),FeatureThres=as.numeric(FeatThres), MitoThres=as.numeric(MtThres))

write.csv(filter.sce, "S02raw_mito_FilterQC.csv")

#Remove Mito genes
#sce <- sce[!is.mito]

#--------------------------Cell Cycle Regression--------------------------

#mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran")) 
#anno <- select(org.Mm.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL") 
#ensembl <- anno$ENSEMBL[match(rownames(sce), anno$SYMBOL)] 
#assignments <- cyclone(sce, mm.pairs, gene.names=ensembl) 

#png("S02_mito_QCplots_CellCycle.png", width =500, height =500, units = "px", pointsize = 17)
#plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)
#dev.off()

#sce <- sce[,assignments$phases=="G1"]

#------------------- Clustering, computing size factor, log normalization ------------------- 

clusters <- quickCluster(sce) 
sce <- computeSumFactors(sce, cluster=clusters)

plot(sizeFactors(sce), sce$nCount_RNA, ylab="Library size", xlab="Size factor")
sce <- logNormCounts(sce, log=FALSE)

#------------------- Convert to Seurat ------------------- 

srt <- as.Seurat(sce)

#------------------- Setting Data w/e sce normalized data ------------------- 

# backup Seurat's norm data
D2 <- NormalizeData(srt)
D2@misc[["seurat_norm_data"]] <- GetAssayData(object = D2, slot = 'data') 

# setting normalized sc2 data
sce.normdata <- as.matrix(log(x = assay(sce, "normcounts") + 1))
sce.normdata [1:10,1:4]

D2@assays[["RNA"]]@data <- sce.normdata

#SetAssayData(object = D2[["RNA"]], slot = "data", new.data = sce.normdata)


D2@assays[["RNA"]]@data[1:10,1:4] #test - is sce data in seurat data slot?
GetAssayData(object = D2, slot = 'data') [1:10,1:4] # same as above

D2@misc[["seurat_norm_data"]] [1:10,1:4] #is seurat normalized data backed up in misc slot


#plotTSNE(sce, colour_by = "Lrrk2")

saveRDS (D2, file = "S02raw_mito_QC.rds")

#---------------------------------------------------------------------------
#-------------------------Legacy------------------------------------------
#---------------------------------------------------------------------------

#Library Size
libsize.drop <- isOutlier(sce$nCount_RNA, nmads=1.0, type="lower", log=TRUE)
LibD <- attr(libsize.drop, "thresholds")
LibD["lower"]
#B02:1.0

#Feature Size
feature.drop <- isOutlier(sce$nFeature_RNA, nmads=0.9, type="lower", log=TRUE)
FeatD <- attr(feature.drop, "thresholds")
FeatD["lower"]
#B02:0.9

# ID cells with high mitochondria gene expression
is.mito <- grepl("^mt-", rownames(sce))
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito)) 

mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=1.5, type="higher")
MitoD <- attr(mito.drop, "thresholds")
MitoD["higher"]
#B02:1.5

#Make QC Plots
png(".png", width =2000, height =500, units = "px", pointsize = 17)
par(mfrow=c(1,4))
#Count
hist(sce$nCount_RNA, xlab="Library size", main="",
     breaks=50, col="grey80", ylab="Number of cells")
abline(v=LibD["lower"], col="blue", lwd=2, lty=2)

hist(sce$nCount_RNA, xlab="Library size", main="",
     breaks=800, col="grey80", xlim=c(0,5000), ylab="Number of cells")
abline(v=LibD["lower"], col="blue", lwd=2, lty=2)
#Features
hist(sce$nFeature_RNA, xlab="Number of expressed genes", main="",
     breaks=100, col="grey80", xlim=c(0,5000), ylab="Number of cells")
abline(v=FeatD["lower"], col="blue", lwd=2, lty=2)
#Mito
hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=40, main="", col="grey80")
abline(v=MitoD["higher"], col="blue", lwd=2, lty=2)
dev.off()

p4 <- plotColData(sce, x="nCount_RNA", y="nFeature_RNA") +
  geom_vline(xintercept=LibD["lower"], color="red", size=1, linetype="dashed") +
  geom_hline(yintercept=FeatD["lower"], color="red", size=1, linetype="dashed")

p5 <- plotColData(sce, x="nCount_RNA", y="pct_counts_Mt") + 
  geom_vline(xintercept=LibD["lower"], color="red", size=1, linetype="dashed")

p6 <- plotColData(sce, x="nFeature_RNA", y="pct_counts_Mt") +
  geom_vline(xintercept=FeatD["lower"], color="red", size=1, linetype="dashed")

png(".png", width =1500, height =500, units = "px", pointsize = 17)
plot_grid(p4,p5,p6,
          nrow=1, ncol=3)
dev.off()

# Remove cells of small lib size, small feature size, high mito content
sce <- sce[,!(libsize.drop | feature.drop)] #add mito.drop if appropriate
filter.sce <- data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), ByMito=sum(mito.drop), Remaining=ncol(sce),
                         CountThres=as.numeric(LibD["lower"]),FeatThres=as.numeric(FeatD["lower"]), MitoThres=as.numeric(MitoD["higher"]))

write.csv(filter.sce, ".csv")
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------