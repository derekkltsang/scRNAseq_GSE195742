BiocManager::install("Seurat")
BiocManager::install("R.utils")
remotes::install_github('satijalab/seurat-wrappers')
remotes::install_github("mojaveazure/seurat-disk")

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(scales)

#SPF_mito_all_traj_v2.rds
#GF_mito_all_traj_v2.rds

#SPF_mito_EC_traj_v2.rds
#GF_mito_EC_traj_v2.rds

#SPF_mito_SEC_traj_v2.rds
#GF_mito_SEC_traj_v2.rds

#GFSPF_mito_ep_labelled_v2.rds

m_srt <- readRDS(file = "GFSPF_mito_ep_labelled_v2.rds")
DimPlot(m_srt, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 1.5)

Idents(m_srt) <- "celltype"
m_srtsub <- subset(m_srt, idents = c("ISC", "TA1", "TA2","V1-V2", 
                                   "V3-V4", "V5-V6"))
saveRDS (m_srt1, file = "GFSPF_mito_EC_traj_v2.rds")

m_srt1 <- subset(m_srt, idents = c("ISC", "TA1", "TA2",
                                   "PC", "GC1", "GC2", 
                                   "EE", "Tuft"))


m_srtsub <- readRDS(file = "GFSPF_mito_SEC_traj_v2.rds")
#GFSPF_mito_EC_traj_v2.rds
#GFSPF_mito_SEC_traj_v2.rds
DimPlot(m_srtsub, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 1.5)

m_srtsub@reductions$umap_original_mito <- m_srtsub@reductions$umap

m_srtsub <- FindVariableFeatures(m_srtsub, selection.method = "vst", nfeatures = 3000)
VariableFeaturePlot(m_srtsub)
m_srtsub <- ScaleData(m_srtsub, verbose = FALSE)
m_srtsub <- RunPCA(m_srtsub, npcs = 30, verbose = FALSE)

m_srtsub <- m_srtsub %>% RunHarmony("orig.ident", plot_convergence = FALSE)
DimHeatmap(m_srtsub, dims = 1:10, cells = 500, balanced = TRUE)
ElbowPlot(m_srtsub)

m_srtsub <- m_srtsub %>% 
  RunUMAP(reduction = "harmony", dims = 1:10) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
  FindClusters(resolution = 0.3, algorithm = 3) %>% 
  identity()

DimPlot(m_srtsub, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 1.5)
saveRDS (m_srtsub, file = "GFSPF_mito_SEC_traj_v2.rds")

#---------------------------------------------------------------------
#build reference datasets
#---------------------------------------------------------------------
Idents(m_srtsub) <- "group"

m_srtsub <- subset(m_srtsub, idents = "SPF")

m_srtsub <- readRDS(file = "SPF_mito_SEC_traj_v2.rds")


test <- substr(colnames(m_srtsub), 1, nchar(colnames(m_srtsub))-4)
m_srtsub <- RenameCells(m_srtsub,new.names = test)

m_srtsub <- RenameCells(m_srtsub, new.names = paste0("S03:",test,"x"))
                     
write.csv(Cells(m_srtsub), file = "SPF_mito_SEC_cellID_obs_v2.csv", row.names = FALSE)
write.csv(Embeddings(m_srtsub, reduction = "umap"), file = "SPF_mito_SEC_cell_embeddings_v2.csv")
write.csv(m_srtsub$celltype, file = "SPF_mito_SEC_clusters_v2.csv")

m_srtsub$celltype

Idents(m_srtsub) <- "celltype"
# Create vector with levels of m_srt@ident
identities <- levels(m_srtsub@active.ident)

# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))

DimPlot(m_srtsub, group.by = "celltype")

colorcsv <- cbind(identities, my_color_palette)

write.csv(colorcsv, file = "SPF_SEC_cluster_colours_v2.csv")



