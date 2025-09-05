#-----------------------------Install required Pkgs-----------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Seurat", dependencies = TRUE)
BiocManager::install("Seurat")

BiocManager::install("dplyr")
BiocManager::install("Matrix")
BiocManager::install("scater")
BiocManager::install("scran")
BiocManager::install("edgeR")
BiocManager::install("ggraph")
BiocManager::install("devtools")
BiocManager::install("clustree")
BiocManager::install("cowplot")
BiocManager::install("Rcpp")
BiocManager::install("monocle")
BiocManager::install("kstreet13/slingshot")
BiocManager::install("topGO")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("plotly")
BiocManager::install("loomR")

BiocManager::install("tidyr")
install.packages("remotes")

remotes::install_github(repo ='mojaveazure/loomR', ref = 'develop')
remotes::install_github("immunogenomics/harmony")


library(clustree)
library(scran)
library(dplyr)
library(scater)
library(Seurat)
library(Matrix)
library(ggplot2)
library(cowplot)
library(loomR)
library(magrittr)
library(harmony)
library(topGO)
library(org.Mm.eg.db)
library(RColorBrewer)
library(plotly)
library(htmlwidgets)
library(viridis)


#------------------- Merge Seurat Objects ------------------- 

srt1 <- readRDS(file = "S03raw_mito_QC.rds")
srt1$orig.ident <- "GF"
srt1$group <- "GF"

srt2 <- readRDS(file = "S02raw_mito_QC.rds")
srt2$orig.ident <- "SPF"
srt2$group <- "SPF"

m_srt <- merge(srt1, srt2, merge.data = TRUE)
rm(srt2)
#------------------- PCA & Integration ------------------- 
DefaultAssay(m_srt) <- "RNA"

m_srt <- FindVariableFeatures(m_srt, selection.method = "vst", nfeatures = 3000)
VariableFeaturePlot(m_srt)
m_srt <- ScaleData(m_srt, verbose = FALSE)
m_srt <- RunPCA(m_srt, npcs = 30, verbose = FALSE)

DimHeatmap(m_srt, dims = 1:25, cells = 500, balanced = TRUE)
#DimPlot(object = m_srt, reduction = "pca", pt.size = .1, group.by = "orig.ident")

m_srt <- m_srt %>% RunHarmony("orig.ident", plot_convergence = FALSE)

#------------------- t-SNE and Clustering ------------------- 
ElbowPlot(m_srt)
m_srt <- m_srt %>% 
  RunUMAP(reduction = "harmony", dims = 1:15) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:15) %>% 
  FindClusters(resolution = 0.3, algorithm = 3) %>% 
  identity()
DimPlot(m_srt, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 1.5)
DimPlot(m_srt, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 1.5)

#------------------- t-SNE and Clustering ------------------- 

#m_srt <- m_srt %>% 
#  RunUMAP(reduction = "pca", dims = 1:20) %>% 
#  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
#  FindClusters(resolution = 0.3, algorithm = 3) %>% 
#  identity()

#------------------- Clustree - optimal resolution ------------------- 

for(x in seq(from=0, to=1, by=0.1)) {
  m_srt <- FindClusters(m_srt, algorithm = 3, resolution = x)}

png("GFSPF_mito_ep_clustree.png", width = 750, height = 750, units = "px", pointsize= 20)
clustree(m_srt)
dev.off()

#select resolution
m_srt <- FindClusters(m_srt, algorithm = 3, resolution = 0.3)

#------------------- Visualization ------------------- 

plot <- DimPlot(m_srt, reduction = "umap")
select.cells <- CellSelector(plot = plot)
m_srt2 <- subset(m_srt, cells = select.cells)

Idents(m_srt) <- "seurat_clusters"

png("GFSPF_mitoQC_ep_UMAP_Nleg_v3.png", width = 800, height = 400, units = "px", res=100)
vp1 <- DimPlot(m_srt, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 1.3) +NoLegend()
vp2 <- DimPlot(m_srt, reduction = "umap", group.by = "orig.ident",label = FALSE, shuffle=TRUE, pt.size = 1.3, cols = c("#F8766D","#078992")) +NoLegend()
plot_grid(vp1, vp2)
dev.off()
?DimPlot
png("GFSPF_mito_ep_UMAP_Nleg_v3.png", width = 650, height = 250, units = "px", res=100)
vp1 <- DimPlot(m_srt, reduction = "umap", group.by = "celltype", label = TRUE,label.size = 3.5, pt.size = 1.5) + theme(text = element_text(size = 12),
                                                                                                                       plot.title = element_blank(),
                                                                                                      axis.title.x=element_blank(),
                                                                                                      axis.text.x=element_blank(),
                                                                                                      axis.title.y=element_blank(),
                                                                                                      axis.text.y=element_blank(),
                                                                                                      plot.margin=unit(c(0.2,2,0.2,0.2),"cm"), 
                                                                                                      legend.key.size = unit(4, 'mm'),
                                                                                                      legend.position=c(1.1, 0.6)) 
vp2 <- DimPlot(m_srt, reduction = "umap", group.by = "orig.ident",label = FALSE, shuffle=TRUE, 
               pt.size = 1.5)+ theme(text = element_text(size = 12),
                                                                    axis.title.x=element_blank(),
                                     plot.title = element_blank(),
                                                                    axis.text.x=element_blank(),
                                                                    axis.title.y=element_blank(),
                                                                    axis.text.y=element_blank(),
                                                                    plot.margin=unit(c(0.2,2,0.2,0.2),"cm"), 
                                                                    legend.position=c(1.0, 0.5))  
plot_grid(vp1, vp2)
dev.off()
#cols = c("#F8766D","#078992")

png("GFSPF_mito_ep_UMAPsplit.png", width = 1000, height = 500, units = "px", res=100)
DimPlot(m_srt, reduction = "umap", split.by = "group", label = TRUE,
        pt.size = 1.3,label.size = 5)
dev.off()

#------------------- Cluster Markers and Cell Type IDs ------------------- 
Idents(m_srt) <- "celltype"
DefaultAssay(m_srt) <- "RNA"
#Find markers for each cluster 
markers <- FindAllMarkers(m_srt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ggtc <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(markers, "GF_SPF_cluster_markers.csv")

DoHeatmap(m_srt, features = c(gtc), group.by = "celltype.group", size = 2, slot = "scale.data", 
          disp.min = -2, disp.max = 2, cells = WhichCells(m_srt, idents = c("6-GC"))) +  
  theme(text = element_text(size = 14))
?DoHeatmap()

#scale_fill_gradientn(colors = c("blue", "white", "red"))
#cells = WhichCells(m_srt, idents = c("5-RSC WT IR", "5-RSC Nod2 IR"))
#top_n(n = 10, wt = -p_val_adj) gives lowest values

genes.to.plot <- c("Lgr5","Ascl2","Axin2","Olfm4", "Slc12a2", "Mki67","Mcm6","Mcm5",
                   "Hes1","Alpi","Reg1" , "Gstm3","Sis", "Slc5a1", "Slc2a2","Apoa4","Ada",
                   "Atoh1","Lyz1","Defa17","Defa22","Defa24","Ang4","Agr2","Muc2","Tff3",
                   "Tph1","Chga","Chgb","Tac1","Scg3","Sct",
                   "Dclk1","Trpm5","Cd3e","Cd3g","Ptprc", "Ms4a1", "mt-Cytb", "mt-Nd1", "mt-Nd3", "mt-Nd4")

genes.to.plot <- c("Lgr5","Ascl2","Axin2","Olfm4", "Slc12a2", "Mki67","Mcm6","Mcm5",
                   "Hes1","Alpi", "Reg1" , "Gstm3","Sis", "Slc5a1", "Slc2a2","Apoa4","Ada",
                   "Atoh1","Lyz1","Defa17","Defa22","Defa24","Ang4","Agr2","Muc2","Tff3",
                   "Tph1","Chga","Chgb","Tac1","Scg3","Sct",
                   "Dclk1","Trpm5")

genes.to.plot <- c("Lgr5","Ascl2","Axin2","Olfm4", "Slc12a2", "Mki67","Mcm6","Mcm5",
                   "Hes1","Alpi", "Reg1" , "Gstm3","Sis", "Slc5a1", "Slc2a2","Apoa4","Ada",
                   "Atoh1","Lyz1","Defa17","Defa22","Defa24","Ang4","Agr2","Muc2","Tff3",
                   "Tph1","Chga","Chgb","Tac1","Scg3","Sct",
                   "Dclk1","Trpm5")
DoHeatmap(m_srt, features = genes.to.plot, group.by = "celltype", size = 2, slot = "scale.data", 
          disp.min = -2, disp.max = 2,)+theme(text = element_text(size = 16))

DoHeatmap(m_srt, features = c(ggtc$gene), group.by = "celltype", size = 2, slot = "scale.data", 
          disp.min = -2, disp.max = 2,)

png("GFSPF_mito_ep_celltype_HM_v2.png", width = 750, height = 750, units = "px", res=100)
DoHeatmap(m_srt, features = genes.to.plot, group.by = "celltype", size = 2, slot = "scale.data", 
          disp.min = -2, disp.max = 2,) +  
  theme(text = element_text(size = 16))
dev.off()
#+scale_fill_viridis()

png("GFSPF_mitoQC_ep_celltype_Dotplot_v2.png", width = 1100, height = 400, units = "px", res=100)
DotPlot(m_srt, features = genes.to.plot) + RotatedAxis() & theme(axis.title.x = element_blank(),axis.title.y = element_blank() )
dev.off()
?DotPlot
#Reorder Groups
#levels(m_srt) <- as.factor(c(0,1,2,3,4,5,6,7,8,9,10))
#m_srt <- m_srt2
#Rename clusters by ID
Idents(m_srt) <- "celltype"
m_srt <- RenameIdents(m_srt, 
                      `0` = "ISC", 
                      `1` = "TA", 
                      `2` = "EPprog",
                      `3` = "V1-2", 
                      `4` = "V3-4",
                      `5` = "V5-6", 
                      `6` = "PC", 
                      `7` = "GC1",
                      `8` = "GC2",
                      `9` = "EE",
                      `10` = "Tuft")

m_srt <- RenameIdents(m_srt, 
                      `ISC` = "ISC", 
                      `TA` = "TA1", 
                      `EPprog` = "TA2",
                      `EC1` = "V1-V2", 
                      `EC2` = "V3-V4",
                      `EC3` = "V5-V6", 
                      `PC` = "PC", 
                      `GC1` = "GC1",
                      `GC2` = "GC2",
                      `EE` = "EE",
                      `Tuft` = "Tuft")

m_srt <- RenameIdents(m_srt, 
                      `ISC` = "ISC", 
                      `TA` = "TA",
                      `EC1` = "V1-V2", 
                      `EC2` = "V3-V4",
                      `EC3` = "V5-V6", 
                      `PC` = "PC", 
                      `GC1` = "GC1",
                      `GC2` = "GC2",
                      `EE` = "EE",
                      `Tuft` = "Tuft",
                      `IM1` = "IM1",
                      `IM2` = "IM2",
                      `Mito` = "Mito")


m_srt$celltype <- Idents(m_srt)
Idents(m_srt) <- factor(Idents(m_srt), 
                            levels = c("ISC", "TA1", "TA2","EC1", 
                                       "EC2", "EC3", 
                                       "PC", "GC1", "GC2", 
                                       "EE", "Tuft"))

Idents(m_srt) <- "celltype.group"
Idents(m_srt) <- factor(Idents(m_srt), 
                        levels = c("ISC_GF","ISC_SPF", "TA_GF","TA_SPF",
                                   "V1-V2_GF", "V1-V2_SPF","V3-V4_GF", "V3-V4_SPF","V5-V6_GF", "V5-V6_SPF", 
                                   "PC_GF", "PC_SPF","GC1_GF", "GC1_SPF","GC2_GF", "GC2_SPF", 
                                   "EE_GF", "EE_SPF","Tuft_GF","Tuft_SPF",  "IM1_GF", "IM1_SPF","IM2_GF","IM2_SPF",
                                   "Mito_GF","Mito_SPF"))

Idents(m_srt) <- "celltype.group"
Idents(m_srt) <- factor(Idents(m_srt), 
                        levels = c("ISC_GF","ISC_SPF", "TA1_GF","TA1_SPF", "TA2_GF", "TA2_SPF",
                                   "V1-V2_GF", "V1-V2_SPF","V3-V4_GF", "V3-V4_SPF","V5-V6_GF", "V5-V6_SPF", 
                                   "PC_GF", "PC_SPF","GC1_GF", "GC1_SPF","GC2_GF", "GC2_SPF", 
                                   "EE_GF", "EE_SPF","Tuft_GF","Tuft_SPF"))
m_srt$celltype.group <- Idents(m_srt)
m_srt$celltype.group <- as.factor(paste(m_srt$celltype, m_srt$group, sep = "_"))
m_srt$celltype <- Idents(m_srt)
cellcount <- table(m_srt$celltype.group)
write.csv(cellcount, file = "GFSPF_mitoQC_ep_cellcount_v2.csv")

# ------------------- make Celltype QC plots  ------------------- 

m_srtqc_sec <- subset(m_srt, idents = c("GC1", "GC2", "PC", "Tuft"))
m_srtqc_ec <- subset(m_srt, idents = c("V1-V2", "V3-V4", "V5-V6"))
m_srtqc_im <- subset(m_srt, idents = c("IM1", "IM2"))
m_srtqc_crypt <- subset(m_srt, idents = c("ISC","TA", "mito"))

png("GFSPF_mito_QCplot_v3.png", width = 1000, height = 1000, units = "px", pointsize= 20)

qcp1 <- FeatureScatter(m_srt, feature1="nCount_RNA", feature2="subsets_Mito_percent", group.by = "celltype") + 
  geom_vline(xintercept=1100, color="red", size=1, linetype="dashed") +
  geom_hline(yintercept=60, color="red", size=1, linetype="dashed") +
  xlim(0,3000) + 
  theme(text = element_text(size = 18), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))

qcp2 <- FeatureScatter(m_srt, feature1="nFeature_RNA", feature2="subsets_Mito_percent", group.by = "celltype") +
  geom_vline(xintercept=600, color="red", size=1, linetype="dashed") +
  geom_hline(yintercept=60, color="red", size=1, linetype="dashed") + 
  theme(text = element_text(size = 18), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))
  
qcp3 <- FeatureScatter(m_srt, feature1="nCount_RNA", feature2="subsets_Mito_percent" , group.by = "group",cols = c("#F8766D","#078992")) +
  geom_vline(xintercept=1100, color="red", size=1, linetype="dashed") +
  geom_hline(yintercept=60, color="red", size=1, linetype="dashed") +
  xlim(0,3000) + 
  theme(text = element_text(size = 18), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))

qcp4 <- FeatureScatter(m_srt, feature1="nFeature_RNA", feature2="subsets_Mito_percent" , group.by = "group",cols = c("#F8766D","#078992")) +
  geom_vline(xintercept=600, color="red", size=1, linetype="dashed") +
  geom_hline(yintercept=50, color="red", size=1, linetype="dashed") + 
  theme(text = element_text(size = 18), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))

plot_grid(qcp1, qcp2, qcp3, qcp4)
dev.off()

saveRDS (m_srt, file = "GFSPF_mito_ep_labelled_v2.rds")

#------------------- Load Dataset ------------------- 
#GFSPF_mito_ep_labelled_v2.rds
m_srt <- readRDS(file = "GFSPF_mito_ep_labelled_v2.rds")
m_srt1 <- readRDS(file = "GFSPF_mitoQC_labelled_epithelium.rds")
m_srt_SC <- readRDS(file = "GF_SPF_SC_v5_filtered.rds")
saveRDS (m_srt, file = "GFSPF_mito_ep_labelled.rds")
#GFSPF_mitoQC_labelled_epithelium.rds
#GF_SPF_labelled_epithelium.rds
#EC_IR_GF_SPF_v2_label.rds
#GF_SPF_labelled.rds
#GF_SPF_IRBS_labelled.rds
#GF_SPF_labelled.rds
#pc_an.rds
#gc_an.rds
#GF_SPF_SC_v5.rds
#------------------- Differential Gene Expression ------------------- 

diff_gene_exp(object = m_srt, cluster_ID = "TA", ident_1 = "SPF", ident_2 = "GF", 
              name = "GFSPF_mito_ep_DE_")

#------------------- Visual Plots ------------------- 
#Violin Plots
DefaultAssay(m_srt) <- "RNA"
Idents(m_srt) <- "celltype"
?VlnPlot
#Mtor, Rptor, Mlst8
png("GFSPF_Vln_Ormdl3.png", width = 450, height = 200, units = "px", pointsize= 20)
VlnPlot(m_srt, features = c("Ormdl3"), group.by = "celltype", split.by = "group",
        pt.size = 0, combine = TRUE) & theme(axis.title.x = element_blank(),
                                             axis.title.y = element_text(size = 20),
                                             axis.text.x = element_text(size = 20), 
                                             axis.text.y = element_text(size = 20), 
                                             text = element_text(size = 23)) 
dev.off()

?VlnPlot()
VlnPlot(m_srtqc_ec, features = c("Reg3b", 'Reg3g', "Dmbt1"), group.by = "celltype", split.by = "group",
        pt.size = 0, combine = TRUE)

VlnPlot(m_srt, features = "Rpl3", group.by = "celltype", split.by = "group",
        pt.size = 0, combine = FALSE)

?VlnPlot

colours.for.plot <- c("#d3d3d3", "#ffeda0" ,"#fed976","#feb24c",
                      "#fd8d3c", "#fc4e2a", "#e31a1c" ,"#bd0026","#800026")

# Gene expression map
png("GFSPF_mitoQC_umap.png", width = 900, height = 250, units = "px", pointsize= 20)
p1 <- FeaturePlot(m_srt, features = c("nFeature_RNA"),
            order = FALSE, min.cutoff = 600, max.cutoff = 1500, 
            cols = c("grey", "red"), pt.size = 1.3) + scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "RdYlBu")))
p2 <- FeaturePlot(m_srt1, features = c("nCount_RNA"),
                  order = FALSE, min.cutoff = 1000, max.cutoff = 5000, 
                  cols = c("grey", "red"), pt.size = 1.3) + scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "RdYlBu")))
p3 <- FeaturePlot(m_srt, features = c("subsets_Mito_percent"),
                  order = FALSE, max.cutoff = 50,
                  cols = c("grey", "red"), pt.size = 1.3) + scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdYlBu")))
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

#rev(brewer.pal(n = 9, name = "RdYlBu"))
#nFeature_RNA
#nCount_RNA
#subsets_Mito_percent

FeaturePlot(m_srt, features = c("Hes1", "Atoh1"),pt.size = 1.0, blend = TRUE, order = TRUE)
?FeaturePlot

FeaturePlot(m_srt, features = c("Alpi"),
            order = FALSE, pt.size = 1.0, split.by = "group") 

# nCount_RNA, nFeature_RNA, pct_counts_Mt

FeaturePlot(m_srt, features = c("mt-Co1"), 
            split.by = "group", 
            max.cutoff =6.0, min.cutoff = 0,
            order = FALSE, pt.size = 1.1)

+ scale_colour_gradientn(colours = colours.for.plot)
P1 + scale_colour_gradientn(colours = colours.for.plot)
#-------------------------------------------------
#---------------- Visualization ------------------ 
#-------------------------------------------------

#m_srt_old <- readRDS(file = "GF_SPF_labelled_v3_overcluster.rds")


Idents(m_srt) <- "group"
m_srt1 <- subset(m_srt, idents = "GF")
m_srt2 <- subset(m_srt, idents = "SPF")

FPlot_split(genetoplot = "Itln1", tf_order = FALSE, maxcutoff = 6.0, name = "GFSPF_mito_ep_")

FeaturePlot(m_srt, features = c("Myd88"), order = FALSE, 
            cols = c("grey", "red"), pt.size = 2) + scale_colour_gradientn(colours = colours.for.plot) + theme(legend.position = "right")

FeaturePlot(m_srt, "Malat1", split.by = "group",
                    min.cutoff = 0, 
                    order = TRUE, pt.size = 1.0) & scale_colour_gradientn(colours = colours.for.plot) & theme(legend.position = "right")

png("GFSPF_mitoEC_Ormdl3_umap.png", width = 500, height = 250, units = "px", pointsize= 20)
FeaturePlot(m_srt, features = c("Ormdl3"), split.by = "group", keep.scale = "feature",
            min.cutoff = 0, max.cutoff = 2.0, order = TRUE, pt.size = 2.0) & 
  scale_colour_gradientn(colours = colours.for.plot) & theme(text = element_text(size = 15), 
                                                             axis.title.x=element_text(size = 15),
                                                             axis.text.x=element_text(size = 15),
                                                             axis.title.y=element_text(size = 15),
                                                             axis.text.y=element_text(size = 15),
                                                             legend.position = c(0.97,0.97),
                                                             legend.text = element_text(size = 14))
dev.off()

png("GFSPF_IR_mitoEC_Ly6a_umap.png", width = 500, height = 250, units = "px", pointsize= 100)
FeaturePlot(m_srt, features = c("Ly6a"), split.by = "group", keep.scale = "feature",
            min.cutoff = 0, max.cutoff = 3.0, order = TRUE, pt.size = 2.0) & 
  scale_colour_gradientn(colours = colours.for.plot) & theme(text = element_text(size = 15),
                                                             axis.title.x=element_blank(),
                                                             axis.text.x=element_blank(),
                                                             axis.title.y=element_blank(),
                                                             axis.text.y=element_blank(),
                                                             legend.position = c(0.97,0.97),
                                                             legend.text = element_text(size = 14))
dev.off()

FeaturePlot(m_srt, features = c("Vegfa"), split.by = "group", keep.scale = "feature",
            min.cutoff = 0,  order = TRUE, pt.size = 2.0) & 
  scale_colour_gradientn(colours = colours.for.plot) & theme(text = element_text(size = 13), 
                                                             legend.position = c(0.95,0.95),
                                                             legend.text = element_text(size = 12),
                                                             axis.title.x=element_blank(),
                                                             axis.text.x=element_blank(),
                                                             axis.title.y=element_blank(),
                                                             axis.text.y=element_blank())

png("GFSPF_mitoEC_umap_Vegfa.png", width = 500, height = 250, units = "px", pointsize= 20)
FeaturePlot(m_srt, features = c("Vegfa"), split.by = "group", keep.scale = "feature",
            min.cutoff = 0, max.cutoff = 2.0, order = TRUE, pt.size = 2.0) & 
  scale_colour_gradientn(colours = colours.for.plot) & theme(text = element_text(size = 15), 
                                                             axis.title.x=element_text(size = 15),
                                                             axis.text.x=element_text(size = 15),
                                                             axis.title.y=element_text(size = 15),
                                                             axis.text.y=element_text(size = 15),
                                                             legend.position = c(0.97,0.97),
                                                             legend.text = element_text(size = 14))
dev.off()

#3D Plot
m_srt_3d <- m_srt

m_srt_3d <- RunUMAP(m_srt_3d, reduction = "harmony", dims = 1:20, n.components = 3L)
Embeddings(object = m_srt_3d, reduction = "umap")

plotting.data <- FetchData(object = m_srt_3d, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "celltype"))

#eye (0-2.5) z = lower view point, x,y = zoom, set closer to origin (x=0, y=0)

#up (0-2.5) >0, set axis as the 'vertical'

#center (set focal point)

scene = list(camera = list(
              eye = list(x = 0.1, y = -1.0, z = 0.15),
             up = list(x = 0.0, y = 0.0, z= 1.0 ),
             center = list(x = 0.0, y = -0.8, z= 0.1)))

plot_ly(data = plotting.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~celltype,
        type = "scatter3d", mode = "markers",
        alpha = 0.8, size = I(100))%>% 
      layout(scene = scene)

plotting.data <- FetchData(object = m_srt, vars = c("UMAP_1", "UMAP_2", "orig.ident"))

fig <- plot_ly(data = plotting.data, x = ~UMAP_1, y = ~UMAP_2, color = ~orig.ident,
               type = "scatter", mode = "markers", colors = c("red","grey"),
               alpha = 0.2, size = I(50), alpha_stroke = 0)
?plot_ly
fig <- fig %>% layout(yaxis = list(zeroline = FALSE, showgrid=FALSE),
                      xaxis = list(zeroline = FALSE, showgrid=FALSE))
fig


#colors = "RdYlBu"
#marker.size = 2

#RdYlBu, YlOrRd

#-------------------------------------------------
#------------------- TopGO ----------------------- 
#-------------------------------------------------

#postive values are increased expression in first group


topGO_diff(cluster_ID = "1", ident_1 = "GF", ident_2 = "SPF", csv_name = "GF_SPF_topGO_1SC_gf.csv")

#-------------------------------------------------
#------------------- Subset ----------------------
#-------------------------------------------------

m_srtsub <- subset(m_srt, idents = c("GF"))
saveRDS (m_srtsub, file = "GF_SPF_SC_v5.rds")
m_srtsub <- readRDS(file = "GFSPF_mito_ISCTAPC_sub.rds")
saveRDS (m_srtsub, file = "GF_mito_Eps.rds")

#-------------------------------------------------
#------------------- Subset ------------------- 
#-------------------------------------------------
DefaultAssay(m_srt) <- "RNA"
Idents(m_srt) <- "celltype"

m_srtsub <- subset(m_srt, idents = c("ISC","TA"))
m_srtsub@reductions$umap_original_mito <- m_srtsub@reductions$umap

plot <- DimPlot(m_srtsub, reduction = "umap_original_mito")
select.cells <- CellSelector(plot = plot)
m_srtsub <- subset(m_srtsub, cells = select.cells)

m_srtsub <- FindVariableFeatures(m_srtsub, selection.method = "vst", nfeatures = 30000)
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
#reduction = "pca"

#------------------------------ Clustree & clustering --------------------------------------
png("GFSPF_mito_ISCTA_clustree.png", width = 750, height = 750, units = "px", pointsize= 20)
for(x in seq(from=0, to=1, by=0.1)) {
  m_srtsub <- FindClusters(m_srtsub, algorithm = 3, resolution = x)}
clustree(m_srtsub)
dev.off()

m_srtsub <- FindClusters(m_srtsub, algorithm = 3, resolution = 0.3)
DimPlot(m_srtsub, reduction = "umap", label = FALSE, pt.size = 1.5, label.size = 6, group.by = "seurat_clusters")
DimPlot(m_srtsub, reduction = "umap", label = FALSE, pt.size = 1.5, label.size = 6, group.by = "orig.ident")

#umap_original_mito
#------------------------------ Plot and labels --------------------------------------------
png("GFSPF_mito_ISCTAPC_UMAP_L_v3.png", width = 600, height = 300, units = "px", res=100)
vp1 <- DimPlot(m_srtsub, reduction = "umap", label = FALSE, pt.size = 1.0, label.size = 6, group.by = "celltype") + NoLegend()
vp2 <- DimPlot(m_srtsub, reduction = "umap", label = FALSE, pt.size = 1.0, label.size = 6, group.by = "group", shuffle=TRUE,cols = c("#F8766D","#078992")) + NoLegend()
plot_grid(vp1, vp2)
dev.off()

genes.to.plot <- c("Lgr5","Olfm4","Ifitm3","Gkn3", "Kcnq1ot1", 
                   "Top2a","Hist1h2ap","Arl6ip1","Mki67","Lbr",
                   "Alpi","Sis", "mt-Cytb", "mt-Nd1","mt-Nd3",
                   "Defa23","Itln1","Lyz1","AY761184","Gm14851")

?DimPlot
Idents(m_srtsub) <-"celltype"
markers <- FindAllMarkers(m_srtsub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
m_srtsub_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

DoHeatmap(m_srtsub, features = genes.to.plot, size=3) + theme(text = element_text(size = 12))

png("GFSPF_mito_ISCTAPC_marker_HM.png", width = 600, height = 300, units = "px", res=100)
DoHeatmap(m_srtsub, features = m_srtsub_markers$gene, size=3) + theme(text = element_text(size = 12))
dev.off()

write.csv(markers, "GFSPF_TA_cluster_markers.csv")

m_srtsub <- RenameIdents(m_srtsub, 
                      `0` = "ISC", 
                      `1` = "TA", 
                      `2` = "EC",
                      `3` = "PC")

Idents(m_srtsub) <- factor(Idents(m_srtsub), 
                        levels = c("ISC", "TA", "EC","PC"))

Idents(m_srtsub) <-"celltype.group"
Idents(m_srtsub) <- factor(Idents(m_srtsub), 
                           levels = c("ISC_GF","ISC_SPF","TA_GF", "TA_SPF", "PC_GF", "PC_SPF","EC_GF", "EC_SPF"))

m_srtsub$celltype.group <- Idents(m_srtsub)

m_srtsub$celltype.group <- as.factor(paste(Idents(m_srtsub), m_srtsub$group, sep = "_"))
m_srtsub$celltype <- Idents(m_srtsub)
cellcount <- table(m_srtsub$celltype.group)
write.csv(cellcount, file = "GFSPF_mito_ISCTAPC_cellcount.csv")

saveRDS (m_srtsub, file = "GFSPF_mito_ISCTAPC_sub.rds")
m_srtsub <- readRDS(file = "GFSPF_mito_ISCTAPC_sub.rds")
#GFSPF_mito_ISC_sub.rds
#GFSPF_mito_PC_sub.rds
#GFSPF_mito_ISCTAcellselect_sub.rds
#GFSPF_mito_TAcellselect_sub.rds
#--------------------------------------------------------------------------------

gtp <- c("Tlr1","Tlr2","Tlr4","Tlr5","Tlr9", 
         "Nod1", "Nod2", "Nlrc4", "Nlrp3", "Nlrp6",
         "Tifa", "Sting", "Tfeb",
         "Myd88", "Ripk2", "Aim2")

gtp <- c("Spink4","Reg3g", "Defa39", "Ang4", 
         "mt-Rnr2", "mt-Co1", "mt-Cytb", "Lyz1",
         "Defa22","Defa29", "Defa30", "Mptx2", 
         "Ssr1",  "Kdelr2", "Defa24")

gtp <- c("Spink4","Reg3g", "Defa39", "Ang4", 
         "mt-Rnr2", "mt-Co1", "mt-Cytb", "Lyz1",
         "Defa22","Defa29", "Defa30", "Mptx2", 
         "Ssr1",  "Kdelr2", "Defa24")

gtp <- c("Agr2", "Reg3g", "Defa39", "Defa24",
         "mt-Rnr2","mt-Co1", "mt-Cytb", "Lyz1",
         "Defa29", "Defa30", "Mptx2")
#"Ang4", "Spink4", "Itln1"

VlnPlot(m_srt, features = gtp, group.by = "celltype", split.by = "group", ncol =4,
        pt.size = 0, combine = TRUE) & theme(axis.title.x = element_blank())

png("GFSPF_ISCTAPC_Mito_Vln.png", width = 450, height = 200, units = "px", pointsize= 20)
VlnPlot(m_srtsub, features = c("subsets_Mito_percent"), group.by = "celltype", split.by = "group",
        pt.size = 0, combine = TRUE) & theme(axis.title.x = element_blank(),
                                             axis.title.y = element_text(size = 20),
                                             axis.text.x = element_text(size = 20), 
                                             axis.text.y = element_text(size = 20), 
                                             text = element_text(size = 23)) + NoLegend()
dev.off()

?VlnPlot
theme(axis.title.x = element_blank())
theme(text = element_text(size = 20))

DoHeatmap(m_srt, features = gtc, group.by = "celltype.group", size = 2, slot = "data", disp.min = 0, disp.max = 1.0) + 
  scale_fill_viridis() +
  theme(text = element_text(size = 20))

F1 <- FeaturePlot(m_srt, features = gtc,
                  min.cutoff = 0, order = TRUE, pt.size = 1.3, combine = FALSE) 

patchwork::wrap_plots(F1) & scale_colour_gradientn(colours = colours.for.plot)

scale_colour_gradientn(colours =rev(brewer.pal(n = 10, name = "RdBu"))) & 
  theme(text = element_text(face = "bold"),
        axis.text.x=element_text(angle=45, hjust=1, size=30),
        axis.title = element_text(size=30,face="bold"),
        axis.title.y.right = element_text(size = 12),
        legend.text=element_text(size=50),
        legend.title=element_text(size=50),
        axis.line = element_line(size=2))


diff_gene_exp(cluster_ID = "4", ident_1 = "4_SPF", ident_2 = "4_GF", 
              name = "GF_SPF_SC_TA_DE_C4")

saveRDS (m_srtsub, file = "GFSPF_SCTA_sub.rds")
m_srtsub <- readRDS(file = "GFSPF_mito_ISCTAPC_sub.rds")
#GF_SPF_SC_v5.rds

#-----------------------------plots----------------------
FeaturePlot(m_srtsub, features = c("Olfm4"), reduction = "umap",
            min.cutoff = 0,
            order = FALSE, pt.size = 2.0) + scale_colour_gradientn(colours = colours.for.plot)
m_srtsub$subsets_Mito_percent

Idents(m_srt) <- "group"
m_srt1 <- subset(m_srt, idents = "GF")
m_srt2 <- subset(m_srt, idents = "SPF")
FPlot_split(genetoplot = "Rps6kb1", tf_order = TRUE, maxcutoff = 2, name = "GFSPF_mito_ep_")

FeaturePlot(m_srt_SC, features = c("Spdef", "Lgr5",
                                     "Hes1"), 
            split.by = "group", 
            max.cutoff = 6, min.cutoff = 0,
            order = FALSE,
            cols = c("grey", "red"), pt.size = 1.5)

#---------------------------------------------------------------------------------------------
#Subset a cluster and relabel with previous labels
m_srt_SC1$sc_idents <- m_srt_SC$seurat_clusters
DimPlot(m_srt_SC, reduction = "umap", group.by = "celltype", label = FALSE, pt.size = 1.5)


plot <- DimPlot(m_srt_SC1, reduction = "umap")
select.cells <- CellSelector(plot = plot)
m_srt_SC1 <- subset(m_srt_SC1, cells = select.cells)

markers <- FindAllMarkers(m_srt_SC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mmtc <- markers %>% group_by(cluster) %>% top_n(n = 12, wt = avg_logFC)

DoHeatmap(m_srt_SC, features = c(mmtc$gene), slot = "data", 
          disp.min = -2, disp.max = 2)

#---------------------------------------------------------------------------------------------
# gene lists

# Dot plot of marker genes

gtc <- c("Ccl5","Gzma","Atf3","Cd74","Ubb","Spink1", "Apoa1","Gip","Fosb","Rgs1",
         "H2-Ab1","Anpep","Ddit3","Klf6","Thbs1","Defa36","Defa23","Mt1","Fabp2")


#C5
gtc <- c("Defa30","Defa24","AY761184","Lyz1","Gzma","Ccl5","Gm14851","Defa29","Apoc3","Lgals4",
         "Fabp2","Rplp1","Defa22","Itln1","Apoa1")

gtc <- c("Defa30","Defa24","AY761184","Lyz1","Gm14851","Defa29","Defa22","Itln1","Gzma","Ccl5","Apoc3","Lgals4",
         "Fabp2","Rplp1","Apoa1", "Atf3", "Fos", "Fosb")

#C6
gtc <- c("Defa30","AY761184","Gm7861","Ccl5","Gzma","Hspa8","Gip","Ubb","Gstp1","Apoa1",
         "Mt2","Cd74","Uba52","Krt19","Krt8","Atf3", "Fos", "Jun")

gtc <- c("Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9", 
                                "Nod1", "Nod2", "Nlrc4", "Nlrp3", "Nlrp6",
                                "Tifa", "Sting","Tfeb",
                                "Myd88", "Ticam1","Ripk2", "Aim2")

gtc <- c("Sectm1a","Sectm1b", "Cd320", "Cd3eap", 
         "Ceacam10", "H2-DMa", "Cd74",
         "H2-Eb1","H2-Aa", "D2-Dmb1", 
         "H2-Ab1",  "Relb")

gtc <- c("Reg3a","Reg3g","Dmbt1","Defa30","Defa24","Defa29","Defa17", "Defa20","Defa3", "Defa23","Itln1", "Lyz1")

F1 <- FeaturePlot(m_srtsub, features = gtc, split.by = "group", keep.scale = "feature",
                  min.cutoff = 0, order = TRUE, pt.size = 2.0) & scale_colour_gradientn(colours = colours.for.plot) & theme(legend.position = "right")

patchwork::wrap_plots(F1) & scale_colour_gradientn(colours = colours.for.plot)

FeaturePlot(m_srtsub, features = c("Nlrc4"), split.by = "group", keep.scale = "feature",
           min.cutoff = 0, max.cutoff = 2.0,order = TRUE, pt.size = 2.5) & 
  scale_colour_gradientn(colours = colours.for.plot) & theme(text = element_text(size = 16), legend.position = c(0.7,0.85),legend.text = element_text(size = 14))

DoHeatmap(m_srt, features = c(gtc), group.by = "celltype", size = 2, slot = "scale.data", 
          disp.min = -2, disp.max = 0.5) +  
  theme(text = element_text(size = 14))

DoHeatmap(m_srt, features = c(gtc), group.by = "celltype", size = 2, slot = "data"
          , disp.max = 3) +  
  theme(text = element_text(size = 14))

VlnPlot(m_srt, features = "Xbp1", group.by = "celltype", split.by = "group",
        pt.size = 0.001, combine = TRUE)& theme(axis.title.x = element_blank())

FeaturePlot(m_srt, features = c("Fzd2"),
            order = TRUE, split.by = "group",
            cols = c("grey", "red"), pt.size = 2) & scale_colour_gradientn(colours = colours.for.plot)& theme(legend.position = "right")

#------------------------
#set order of clusters
m_srt <- RenameIdents(m_srt, 
                      `0` = "0-SC", 
                      `1` = "1-EC", 
                      `2` = "2-EPprog",
                      `3` = "3-TA", 
                      `4` = "4-EC", 
                      `5` = "5-GC",
                      `6` = "6-EC", 
                      `7` = "7-GC", 
                      `8` = "8-PC", 
                      `9` = "9-EE",
                      `10` = "10-Tuft")

Idents(m_srt) <- factor(Idents(m_srt), 
                        levels = c("0-SC", "1-TA", "2-EPprog", 
                                   "3-EC", "4-EC", "5-EC", 
                                   "6-EC", "7-PC", "8-PC", 
                                   "9-EE", "10-TC"))
