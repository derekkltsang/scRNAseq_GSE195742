source("http://bioconductor.org/biocLite.R")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
#install.packages("devtools")
BiocManager::install("terra")
devtools::install_github('cole-trapnell-lab/leidenbase', force = TRUE)
devtools::install_github('cole-trapnell-lab/monocle3')

library(monocle3)

m_srt <- readRDS(file = "GFSPF_mito_ep_labelled_v2.rds")
Idents(m_srt) <- "group"
m_srt <- subset(m_srt, idents = c("GF"))

Idents(m_srt) <- "celltype"
m_srt1 <- subset(m_srt, idents = c("ISC", "TA1", "TA2","V1-V2", 
                                   "V3-V4", "V5-V6"))
saveRDS (m_srt1, file = "GF_mito_EC_traj_v2.rds")

m_srt1 <- subset(m_srt, idents = c("ISC", "TA1", "TA2",
                                   "PC", "GC1", "GC2", 
                                   "EE", "Tuft"))
saveRDS (m_srt1, file = "GF_mito_SEC_traj_v2.rds")
#SPF_mito_all_traj_v2.rds
#GF_mito_all_traj_v2.rds

#SPF_mito_EC_traj_v2.rds
#GF_mito_EC_traj_v2.rds

#SPF_mito_SEC_traj_v2.rds
#GF_mito_SEC_traj_v2.rds
m_srt <- readRDS(file = "GF_mito_EC_traj_v2.rds")

#Extract count data, phenotype data, and feature data from the Seurat Object.
counts.data <- as(as.matrix(m_srt@assays$RNA@data), 'sparseMatrix')
#counts.data <- as(as.matrix(m_srt@assays$RNA@counts), 'sparseMatrix')
pheno.data <- m_srt@meta.data
feature.data <- data.frame(gene_short_name = row.names(counts.data), row.names = row.names(counts.data))

#Construct a CellDataSet.
cds <- new_cell_data_set(expression_data = counts.data, 
                         cell_metadata = pheno.data, 
                         gene_metadata = feature.data)


cds <- preprocess_cds(cds, num_dim = 4, norm_method = "none",pseudo_count=0)
plot_pc_variance_explained(cds)
?preprocess_cds
#cds <- preprocess_cds(cds, num_dim = 100, norm_method = "log")

#?preprocess_cds
#plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds, max_components = 2,
                        reduction_method = "UMAP",
                        preprocess_method = "PCA")
#?reduce_dimension
plot_cells(cds, color_cells_by="celltype",   cell_size = 1.0,   group_label_size = 5, reduction_method = "UMAP")
#cds <- cluster_cells(cds, resolution=1e-2)
cds <- cluster_cells(cds, resolution=0.0002, reduction_method = "UMAP")

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by="celltype", cell_size = 1.0, group_label_size = 5, group_cells_by = c("cluster"))
plot_cells(cds, color_cells_by="celltype", label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 1.0, group_label_size = 5)
