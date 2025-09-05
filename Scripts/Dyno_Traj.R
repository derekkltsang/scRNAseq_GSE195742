#install.packages("devtools")
devtools::install_github("dynverse/dyno")
BiocManager::install("tidyverse")
BiocManager::install("fontawesome")
BiocManager::install("devtools")
BiocManager::install("dynwrap")
install.packages("hdf5r")

dynwrap::test_docker_installation(detailed = TRUE)

#all cells colours
your_colors <- c('#EF67EB','#00A6FF', '#B385FF',
                 '#F8766D', '#00BADE', '#DB8E00', 
                 '#AEA200', '#FF63B6', '#64B200',
                 '#00BD5C','#00C1A7')

#SEC only
your_colors <- c('#EF67EB','#00A6FF', '#B385FF',
                 '#F8766D', '#00BADE', '#DB8E00', 
                 '#AEA200', '#FF63B6')


library(dyno)
library(tidyverse)
library(fontawesome)
library(devtools)
library(dynwrap)
#-----------------------------------------------------------------------------------------------------
m_srt <- readRDS(file = "GFSPF_mito_ep_labelled_v2.rds")
Idents(m_srt) <- "group"
m_srt <- subset(m_srt, idents = "SPF")
#-----------------------------------------------------------------------------------------------------
m_srt <- readRDS(file = "GFSPF_mito_ep_labelled_v2.rds")
DimPlot(m_srt, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 1.5)

Idents(m_srt) <- "celltype"
m_srtsub <- subset(m_srt, idents = c("ISC", "TA1", "TA2","V1-V2", 
                                     "V3-V4", "V5-V6"))
#saveRDS (m_srtsub, file = "GFSPF_mito_EC_traj_v2.rds")

m_srtsub <- subset(m_srt, idents = c("ISC", "TA1", "TA2",
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

Idents(m_srtsub) <- "group"
m_srtsub1 <- subset(m_srtsub, idents = "SPF")
DimPlot(m_srtsub1, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 1.5)

saveRDS (m_srtsub1, file = "SPF_mito_SEC_traj_v2.rds")
########################################################### SPF #############################################################

#SPF_mito_all_traj_v2.rds
#GF_mito_all_traj_v2.rds

#SPF_mito_EC_traj_v2.rds
#GF_mito_EC_traj_v2.rds

#SPF_mito_SEC_traj_v2.rds
#GF_mito_SEC_traj_v2.rds

#GFSPF_mito_ep_labelled_v2.rds
srt_spf <- readRDS(file = "SPF_mito_SEC_traj_v2.rds")
DimPlot(srt_spf, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 1.5) +scale_color_manual(values=your_colors)


object_counts_spf <- Matrix::t(as(as.matrix(srt_spf@assays$RNA@counts), 'sparseMatrix'))
object_expression_spf <- Matrix::t(as(as.matrix(srt_spf@assays$RNA@data), 'sparseMatrix'))

object_dimred_spf <- as.matrix(srt_spf@reductions$umap@cell.embeddings)

object_dyn_spf<- wrap_expression(
  counts = object_counts_spf, 
  expression = object_expression_spf)

object_dyn_spf <- add_grouping(
  object_dyn_spf, srt_spf$celltype)

object_dyn_spf <- add_dimred(
  object_dyn_spf, object_dimred_spf)
rm(start_cells)
Idents(srt_spf) <- "celltype"
start_cells_spf <- WhichCells(srt_spf, idents= "ISC")
object_dyn_spf <- add_prior_information(
  object_dyn_spf, start_id = start_cells_spf)

#guidelines_shiny(dataset = object_dyn_spf)

#Dyno model selects which of the >50 models are best for your data
methods_selected <- guidelines$methods_selected

#Build model
#model_ss <- infer_trajectory(object_dyn_spf, method = "slingshot")
#model_paga <- infer_trajectory(object_dyn_spf, method = "paga_tree")

model_paga_spf <- infer_trajectory(object_dyn_spf, ti_paga_tree(n_comps = as.integer(5)))
#GF SEC - as.integer(10) SPF SEC - as.integer(5)
#GF all - as.integer(150) SPF all - as.integer(150) 
#model <- model %>% add_dimred(dyndimred::dimred_mds, expression_source = object_dyn_spf$expression)
model_paga_spf <- model_paga_spf %>% add_dimred(object_dimred_spf, expression_source = object_dyn_spf$expression)
model_paga_rooted_spf <- model_paga_spf %>% add_root_using_expression("Lgr5", expression_source = object_dyn_spf)

png("SPF_mito_SEC_pagatree.png", width = 450, height = 450, units = "px", pointsize= 20)
plot_dimred(
  model_paga_rooted_spf, size_milestones = 3, size_transitions =1.5, size_cells	= 2.5,
  expression_source = object_dyn_spf$expression,   arrow = grid::arrow(type = "closed", length = unit(0.15, "inches")),
  grouping = srt_spf$celltype) + scale_color_manual(values=your_colors) + theme(legend.position = "none")
dev.off()
?plot_dimred()
png("SPF_mito_SEC_dendo.png", width = 800, height = 350, units = "px", pointsize= 100)
plot_dendro(model_paga_rooted_spf, size_cells = 3.5,
            grouping = srt_spf$celltype)+ scale_color_manual(values=your_colors) + theme(legend.position = "none")
dev.off()

plot_dimred(model_paga_rooted_spf, feature_oi = "Sox4", expression_source = object_dyn_spf) + ggtitle("Feature expression")

########################################################### GF #############################################################

srt_gf <- readRDS(file = "GF_mito_SEC_traj_v2.rds")
DimPlot(srt_gf, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 1.5)


object_counts_gf <- Matrix::t(as(as.matrix(srt_gf@assays$RNA@counts), 'sparseMatrix'))
object_expression_gf <- Matrix::t(as(as.matrix(srt_gf@assays$RNA@data), 'sparseMatrix'))

object_dimred_gf <- as.matrix(srt_gf@reductions$umap@cell.embeddings)

object_dyn_gf<- wrap_expression(
  counts = object_counts_gf, 
  expression = object_expression_gf)

object_dyn_gf <- add_grouping(
  object_dyn_gf, srt_gf$celltype)

object_dyn_gf <- add_dimred(
  object_dyn_gf, object_dimred_gf)
rm(start_cells)
Idents(srt_gf) <- "celltype"
start_cells_gf <- WhichCells(srt_gf, idents= "ISC")
object_dyn_gf <- add_prior_information(
  object_dyn_gf, start_id = start_cells_gf)

#guidelines_shiny(dataset = object_dyn_gf)

#Dyno model selects which of the >50 models are best for your data
methods_selected <- guidelines$methods_selected

#Build model
#model_ss <- infer_trajectory(object_dyn_gf, method = "slingshot")
#model_paga <- infer_trajectory(object_dyn_gf, method = "paga_tree")

model_paga_gf <- infer_trajectory(object_dyn_gf, ti_paga_tree(n_comps = as.integer(10)))
#GF SEC - as.integer(10) gf SEC - as.integer(5)
#GF all - as.integer(150) gf all - as.integer(150) 50
#model <- model %>% add_dimred(dyndimred::dimred_mds, expression_source = object_dyn_gf$expression)
model_paga_gf <- model_paga_gf %>% add_dimred(object_dimred_gf, expression_source = object_dyn_gf$expression)
model_paga_rooted_gf <- model_paga_gf %>% add_root_using_expression("Lgr5", expression_source = object_dyn_gf)

png("GF_mito_SEC_pagatree.png", width = 450, height = 450, units = "px", pointsize= 20)
plot_dimred(
  model_paga_rooted_gf, size_milestones = 3, size_transitions =1.5, size_cells	= 2.5,
  expression_source = object_dyn_gf$expression,   arrow = grid::arrow(type = "closed", length = unit(0.15, "inches")),
  grouping = srt_gf$celltype) + scale_color_manual(values=your_colors) + theme(legend.position = "none")
dev.off()

png("GF_mito_SEC_dendo.png", width = 800, height = 350, units = "px", pointsize= 100)
plot_dendro(model_paga_rooted_gf,   size_cells = 3.5,
            grouping = srt_gf$celltype)+ scale_color_manual(values=your_colors)+ theme(legend.position = "none")
dev.off()

?plot_dendro
plot_dimred(model_paga_rooted_gf, feature_oi = "Sox4", expression_source = object_dyn_gf) + ggtitle("Feature expression")
###################################################################################################################################

plot_dimred(
  model_paga_rooted_spf, size_milestones = 3, size_transitions =1, size_cells	= 2.5,
  expression_source = object_dyn_spf$expression, 
  grouping = srt_spf$celltype) + scale_color_manual(values=your_colors)


plot_dimred(model_paga_rooted_spf, feature_oi = "Sox4", expression_source = object_dyn_spf) + ggtitle("Feature expression")

plot_dendro(model_paga_rooted_spf, 
            grouping = srt_spf$celltype)

#all cells colours
your_colors <- c('#EF67EB','#00A6FF', '#B385FF',
                 '#F8766D', '#00BADE', '#DB8E00', 
                 '#AEA200', '#FF63B6', '#64B200',
                 '#00BD5C','#00C1A7')

#SEC only
your_colors <- c('#EF67EB','#00A6FF', '#B385FF',
                 '#F8766D', '#00BADE', '#DB8E00', 
                 '#AEA200', '#FF63B6')
your_colors <- c('#F8766D','#DB8E00', '#AEA200',
                 '#00BADE', '#00A6FF', '#B385FF', 
                 '#EF67EB', '#FF63B6')

###########################################################################################

ISC	#F8766D
TA1	#DB8E00
TA2	#AEA200
V1-V2	#64B200
V3-V4	#00BD5C
V5-V6	#00C1A7
PC	#00BADE
GC1	#00A6FF
GC2	#B385FF
EE	#EF67EB
Tuft	#FF63B6


# Reproduces the guidelines as created in the shiny app - SPF ALL
answers <- dynguidelines::answer_questions(
  multiple_disconnected = NULL, 
  expect_topology = NULL, 
  expected_topology = NULL, 
  n_cells = 3272, 
  n_features = 15386, 
  memory = "30GB", 
  prior_information = "start_id", 
  docker = TRUE
)
guidelines <- dynguidelines::guidelines(answers = answers)

# Reproduces the guidelines as created in the shiny app - GF ALL
answers <- dynguidelines::answer_questions(
  multiple_disconnected = NULL, 
  expect_topology = NULL, 
  expected_topology = NULL, 
  n_cells = 2365, 
  n_features = 15386, 
  memory = "30GB", 
  prior_information = "start_id", 
  docker = TRUE
)
guidelines <- dynguidelines::guidelines(answers = answers)