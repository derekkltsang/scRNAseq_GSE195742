if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("AnnotationHub")
BiocManager::install("ComplexHeatmap")


library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationHub)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(ComplexHeatmap)

library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(digest)
library(cluster)
library(viridis)

#GFSPF_mito_ep_labelled_v2
#EC_IR_GF_SPF_v2_label.rds
#GFSPF_IR_EC_SCT.rds
m_srt <- readRDS(file = "GFSPF_mito_ep_labelled_v2.rds")
Idents(m_srt) <- "celltype"

DimPlot(m_srt, reduction = "umap", group.by = "celltype",label = TRUE, label.size = 6)

m_srt.markers <- FindAllMarkers(m_srt, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
m_srt.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

##################################################################
#Subsetting top 100 markers with adjusted p values lower than .05#
##################################################################
top100 <- m_srt.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)

df <- top100pval[,7:6]
dfsample <- split(df$gene,df$cluster)
length(dfsample)

#The output of length(dfsample) returns how many clusters you have
#Here there at 9 clusters (0, 1, 2, 3, 4, 5, 6, 7 and 8)
#I'm sure there's a better way but you have to make a line like below for each cluster

dfsample$`ISC` = bitr(dfsample$`ISC`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`TA` = bitr(dfsample$`TA`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`EPprog` = bitr(dfsample$`EPprog`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`EC1` = bitr(dfsample$`EC1`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`EC2` = bitr(dfsample$`EC2`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`EC3` = bitr(dfsample$`EC3`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`PC` = bitr(dfsample$`PC`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`GC1` = bitr(dfsample$`GC1`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`GC2` = bitr(dfsample$`GC2`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`EE` = bitr(dfsample$EE, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`Tuft` = bitr(dfsample$`Tuft`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

#do the same here, a line like below for each cluster
genelist <- list("ISC" = dfsample$`ISC`$ENTREZID, 
                 "TA" = dfsample$`TA`$ENTREZID,
                 "EPprog" = dfsample$`EPprog`$ENTREZID,
                 "EC1" = dfsample$`EC1`$ENTREZID,
                 "EC2" = dfsample$`EC2`$ENTREZID,
                 "EC3" = dfsample$`EC3`$ENTREZID,
                 "PC" = dfsample$`PC`$ENTREZID,
                 "GC1" = dfsample$`GC1`$ENTREZID,
                 "GC2" = dfsample$`GC2`$ENTREZID,
                 "EE" = dfsample$`EE`$ENTREZID,
                 "Tuft" = dfsample$`Tuft`$ENTREZID)

GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db")
write.csv(GOclusterplot, file = "GFSPF_Celltype_GOcluster.csv")
dotplot(GOclusterplot)
cnetplot(GOclusterplot)

enrichWPclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichWP", organism = "Mus musculus")
write.csv(enrichWPclusterplot, file = "GFSPF_Celltype_WPcluster.csv")
dotplot(enrichWPclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG")
dotplot(KEGGclusterplot)

##################################################################
#GO enrichment analysis
##################################################################
Idents(m_srt) <- "celltype"
DimPlot(m_srt, reduction = "umap", group.by = "celltype",label = TRUE, label.size = 6)

m_srt1 <- subset(m_srt, idents = c("PC"))
Idents(m_srt1) <- "group"

m_srt1.markers <- FindMarkers(m_srt1, ident.1 = "SPF", ident.2 = "GF", only.pos = FALSE, verbose = FALSE)
m_srt1.markers  %>% top_n(n = 2, wt = avg_log2FC)
m_srt1.markers$cluster <- ifelse(m_srt1.markers$avg_log2FC>0, "SPF", "GF")

write.csv(m_srt1.markers, file = "GFSPF_mito_ep_DE_PC.csv")

pval <- subset(m_srt1.markers, rowSums(m_srt1.markers[5] < 0.05) > 0)
pval$cluster <- ifelse(pval$avg_log2FC>0, "SPF", "GF")
pval$gene <- rownames(pval)

df1 <- pval[,7:6]
df1sample <- split(df1$gene,df1$cluster)
length(df1sample)

df1sample$`SPF` = bitr(df1sample$`SPF`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
df1sample$`GF` = bitr(df1sample$`GF`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

genelist1 <- list("SPF" = df1sample$`SPF`$ENTREZID, 
                 "GF" = df1sample$`GF`$ENTREZID)

GOclusterplot <- compareCluster(geneCluster = genelist1, fun = "enrichGO", OrgDb = "org.Mm.eg.db")
write.csv(GOclusterplot, file = "GFSPF_mito_EC_GOcluster.csv")
?dotplot
dotplot(GOclusterplot)
cnetplot(GOclusterplot)

enrichWPclusterplot <- compareCluster(geneCluster = genelist1, fun = "enrichWP", organism = "Mus musculus")
write.csv(enrichWPclusterplot, file = "GFSPF_mito_EC_WPcluster.csv")
dotplot(enrichWPclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = genelist1, fun = "enrichKEGG")
dotplot(KEGGclusterplot)

######################################################################
#volcano plots
######################################################################
m_srt1 <- subset(m_srt, idents = c("EC1","EC2","EC3"))
Idents(m_srt) <- "group"
m_srt1.markers <- FindMarkers(m_srt, ident.1 = "SPF", ident.2 = "GF", only.pos = FALSE,
                              logfc.threshold=0, min.pct=0, min.cells.group = 1, min.cells.feature = 1)
m_srt1.markers$cluster <- ifelse(m_srt1.markers$avg_log2FC>0, "SPF", "GF")
m_srt1.markers$gene <- rownames(m_srt1.markers)

de <- as.data.frame(m_srt1.markers[,c("gene","avg_log2FC","p_val_adj", "cluster")]) 

#log2FoldChange.SPF_SIvsGF_SI.
#pvalue.SPF_SIvsGF_SI.
#padj.SPF_SIvsGF_SI.
#significant.SPF_SIvsGF_SI.

#labels based on sig from fold exp and pvalue
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$avg_log2FC > 0.75 & de$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$avg_log2FC < -0.75 & de$p_val_adj < 0.05] <- "DOWN"

de$delabel[de$diffexpressed != "NO"] <- de$gene[de$diffexpressed != "NO"]

# plots
p1 <- ggplot(de, aes(x=avg_log2FC, y=-log(p_val_adj))) + geom_point() + theme_minimal()

p2 <- p1 + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

ggplot(de, aes(x=avg_log2FC, y=-log(p_val_adj), label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = Inf,label.padding= 0.5, force=5)

?geom_text_repel
##############################################
png("GFSPF_bulkRNA_ISC_DEs_small.png", width = 400, height = 400, units = "px", res=100)
ggplot(de, aes(x=avg_log2FC, y=-log10(p_val_adj),
               col=delabel,
               label=delabel)) + 
  geom_point(size = 2) + 
  geom_text_repel(max.overlaps = Inf, direction = "y", nudge_x = 2) +
  scale_color_manual(values=c("grey", "black")) +
  theme_minimal() +
  xlab("Log 2 Fold Change")+
  ylab("-log10(adj pvalue)")+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        legend.position = "none") +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
dev.off()



######################################################################
#volcano plots
######################################################################
Idents(m_srt) <- "celltype"

m_srt1 <- subset(m_srt, idents = c("PC", "ISC", "TA"))
m_srt1 <- m_srt
Idents(m_srt1) <- "group"

m_srt1.markers <- FindMarkers(m_srt1, ident.1 = "SPF", ident.2 = "GF", only.pos = FALSE,
                              logfc.threshold=0, min.pct=0, min.cells.group = 1, min.cells.feature = 1)
m_srt1.markers$cluster <- ifelse(m_srt1.markers$avg_log2FC>0, "SPF", "GF")
m_srt1.markers$gene <- rownames(m_srt1.markers)
write.csv(m_srt1.markers, file = "GFSPF_mito_ep_DEallgenes_Crypt.csv")

de <- as.data.frame(m_srt1.markers[,c("gene","avg_log2FC","p_val_adj", "cluster")]) 
#labels based on sig from fold exp and pvalue
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$avg_log2FC > 0.75 & de$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$avg_log2FC < -0.75 & de$p_val_adj < 0.05] <- "DOWN"

de$delabel[de$diffexpressed != "NO"] <- de$gene[de$diffexpressed != "NO"]
gtc <- de$delabel

png("GFSPF_mito_AllEps_StarvTFs_Volcano.png", width = 500, height = 600, units = "px", res=100)
EnhancedVolcano(de,
                lab = de$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = c(gtc),
                drawConnectors = TRUE,
                FCcutoff = 0.75,
                pCutoff = 0.05,
                labSize = 4,
                maxoverlapsConnectors = Inf)
dev.off()
#xlim = c(-log2(10e2), -log2(10e-2)),
?EnhancedVolcano

#-------------------------------------

gtc <- de$delabel
gtc <- c("Reg3a","Reg3g","Dmbt1","Defa30","Defa24","Defa29","Defa17", "Defa20","Defa3", "Defa23","Itln1", "Lyz1")
gtc <- c("Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9", 
         "Nod1", "Nod2", "Nlrc4", "Nlrp3", "Nlrp6",
         "Tifa", "Sting","Tfeb",
         "Myd88", "Ticam1","Ripk2", "Aim2")

gtc <- c("Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9", 
         "Nod1", "Nod2", "Nlrc4", "Nlrp3", "Nlrp6",
         "Tifa", "Sting","Tfeb",
         "Myd88", "Ticam1","Ripk2", "Aim2")

gtc <- c("Tlr1","Tlr2","Tlr3","Tlr6","Tlr7", 
         "Nod1", "Nod2", "Nlrc4", "Nlrp6",
         "Myd88", "Ticam1","Ripk2", "Aim2")

gtp <- c("Tgfbr2","Tgfbr1","Il22ra1","Il6ra","Il4ra","Il1rap","Il1r1","Cxcr2","Ifnlr1","Ifnar2","Ifnar1","Il18r1", 
         "Tnfrsf1a", "Ifngr2","Ifngr1")
gtp <- c("Aim2","Ripk2","Ticam1","Myd88","Nlrp6","Nlrc4", "Nod2",
         "Nod1","Tlr4","Tlr3","Tlr2","Tlr1")

gtc <- c("Aim2","Ripk2","Ticam1","Myd88","Nlrp6","Nlrc4", "Nod2",
         "Nod1","Tlr5","Tlr4","Tlr3","Tlr2","Tlr1")

gtc <- c("Idh1","Ldha","Adh1","Gpd1","Akr1b7","Hpgd","Sord","Akr1a1","Mdh2", 
         "Cbr1", "Ehhadh", "Mdh1", "Hadh", "Ugdh",
         "Akr1c13", "Me2","Akr1c12",
         "Cryl1", "Adh6a","Akr1c19")

gtc <- c("Atf3", "Atf5", "Egr2", "Egr3", "Egr4", "Fos", "Fosb", "Fosl1", "Foxq1", 
       "Gbx1", "Jun", "Prdx1", "Klf6", "Maff", "Nr4a3", "Rax", "Srf", "Tceal7", 
       "Vgll3", "Zmynd15", "Ccl20", "Cyr61", "Cxcl1", "Cxcl10", "Edn1", "Edn2", "Epha7", 
       "F3", "Flrt3", "Hbegf", "Hgf", "Itga1", "Jam3", "Plau", "Plaur", "Thbs1")

gtc <- c("Atf3", "Atf5", "Egr2", "Egr3", "Egr4", "Fos", "Fosb", "Fosl1", "Foxq1", 
         "Gbx1", "Jun", "Prdx1", "Klf6", "Maff", "Nr4a3", "Rax", "Srf", "Tceal7", 
         "Vgll3", "Zmynd15")

gtc <- c("Atf3", "Atf5", "Egr2", "Egr4", "Fos", "Fosb", "Fosl1", "Foxq1", 
         "Gbx1", "Jun", "Prdx1", "Klf6", "Maff", "Nr4a3", "Rax", "Srf", "Tceal7", 
         "Vgll3", "Ccl20", "Cyr61", "Cxcl1", "Edn1", "Edn2", "Hbegf", "Hgf", "Plau", "Thbs1")

gtc <- c("Ifit1","Ifit3","Xaf1","Stat1","Oas1a","Bst2","Igtp","Cdc34","Irgm1","Ifitm3","Oas1g","Ifi47")

?EnhancedVolcano

EnhancedVolcano(de,
                lab = de$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = c(gtc),
                drawConnectors = TRUE,
                FCcutoff = 0.75,
                pCutoff = 0.05,
                labSize = 4,
                maxoverlapsConnectors = Inf)

FeaturePlot(m_srt, features = c("Tlr1"),
            order = TRUE, split.by = "group",   max.cutoff = 7,
            cols = c("grey", "red"), pt.size = 2) & scale_colour_gradientn(colours = colours.for.plot)& theme(legend.position = "right")

FeaturePlot(m_srt, features = c("Fabp2"),
            order = FALSE, split.by = "group",
            cols = c("grey", "red"), pt.size = 2) & scale_colour_gradientn(colours = colours.for.plot)& theme(legend.position = "right")

VlnPlot(m_srt, features = c("Rps6kb1", "Rps6"), group.by = "celltype", split.by = "group",
        pt.size = 0, combine = TRUE)
###########################################  Proportion of 'gene' expressing cells in cluster ########################################### 

Idents(m_srt) <- "celltype"
#features must be present in at least one cell or else DE crashes
m_srt1.pct <- FindAllMarkers(m_srt, features = gtp,
                           logfc.threshold=0, min.pct=0, min.cells.group = 1, min.cells.feature = 1)

m_srt1.pct <- as.matrix(PrctCellExpringGene(object = m_srt, genes = gtc, group.by = "celltype.group"))

#plot from tidy data
m_srt1.pct <- as.data.frame(PrctCellExpringGene(object = m_srt, genes = gtp, group.by = "celltype.group"))
gtc <- c("Aim2","Ripk2","Ticam1","Myd88","Nlrp6","Nlrc4", "Nod2",
         "Nod1","Tlr4","Tlr3","Tlr2","Tlr1")

gtc <- c("Tlr1","Tlr2","Tlr3","Tlr4","Myd88","Ticam1")
gtc <- c("Nod1","Nod2","Ripk2")
gtc <- c("Nlrp6","Nlrc4","Pycard","Casp1")
m_srt1.pct <- as.data.frame(PrctCellExpringGene(object = m_srt, genes = gtc, group.by = "celltype"))

m_srt1.pct$Markers <- factor(m_srt1.pct$Markers, levels=gtc)

m_srt1.pct$Feature <- factor(m_srt1.pct$Feature, levels=c("RevSC","ISC","EC1","EC2","EC3","GCPC","EE"))
m_srt1.pct$Feature <- factor(m_srt1.pct$Feature, levels=c("RevSC_GF","RevSC_SPF","ISC_GF","ISC_SPF","EC1_GF","EC1_SPF",
                                                          "EC2_GF","EC2_SPF","EC3_GF","EC3_SPF","GCPC_GF","GCPC_SPF","EE_GF","EE_SPF"))

m_srt1.pct$Feature <- factor(m_srt1.pct$Feature, levels=c("RevSC","ISC","TA","EC1","EC2","EC3","EC4","EC5","GC","PC","EEC1","EEC2"))
m_srt1.pct$Feature <- factor(m_srt1.pct$Feature, levels=c("ISC_GF","ISC_SPF", "TA1_GF","TA1_SPF", "TA2_GF", "TA2_SPF",
                                                          "V1-V2_GF", "V1-V2_SPF","V3-V4_GF", "V3-V4_SPF","V5-V6_GF", "V5-V6_SPF", 
                                                          "PC_GF", "PC_SPF","GC1_GF", "GC1_SPF","GC2_GF", "GC2_SPF", 
                                                          "EE_GF", "EE_SPF","Tuft_GF","Tuft_SPF"))
#("RevSC","ISC","TA","EC1","EC2","EC3","EC4","EC5","GC","PC","EEC1","EEC2")
m_srt1.pct$Feature <- factor(m_srt1.pct$Feature, levels=c("ISC_GF","ISC_SPF","TA_GF", "TA_SPF", "PC_GF", "PC_SPF","EC_GF", "EC_SPF"))

m_srt1.pct$Feature <- factor(m_srt1.pct$Feature, levels=c("RevSC1","RevSC2", "ISC","EC1","EC2", "EC3","GC", "PC","EE","Tuft"))

m_srt1.pct$Feature <- factor(m_srt1.pct$Feature, levels=c("GF CTRL","SPF CTRL","GF IR", "SPF IR"))
#m_srt1.pct$Cell_proportion <- factor(m_srt1.pct$Cell_proportion)

#DoHeatmap(m_srt, features = gtp, group.by = "celltype.group", size=2, slot= "data") + theme(text = element_text(size = 13))


png("GFSPF_IRBS_eps_NODs_HM_CellProp_lim10_v3.png", width = 700, height = 310, units = "px", res=100)
ggplot(m_srt1.pct, aes(Feature, Markers)) + geom_tile(aes(fill=Cell_proportion)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0, size =20),
        axis.text.y = element_text(size =20)) + 
  scale_x_discrete(position = "top") + scale_fill_viridis(limits = c(0, 10),na.value = "yellow")+ 
  theme(text = element_text(size = 20))
dev.off()

#geom_tile(aes(fill=Cell_proportion),colour="white",size=0.25)
#scale_fill_manual(values=c("#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598","#abdda4","#ddf1da"))
#scale_fill_viridis(limits=c(0,0.2))

######################################## #complexHM ######################################## #
m_srt1.pct <- as.data.frame(PrctCellExpringGene(object = m_srt1, genes = gtc, group.by = "celltype.group"))
m_srt1.pctgenes <- m_srt1.pct$Markers
m_srt1.pct <- pivot_wider(m_srt1.pct, names_from = Feature, values_from = Cell_proportion)
m_srt1.pct <- data.matrix(m_srt1.pct)

rownames(m_srt1.pct) <- m_srt1.pctgenes
m_srt1.pct<- m_srt1.pct[, colnames(m_srt1.pct) != "Markers"]


myCol <- colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
myBreaks <- seq(-3, 3, length.out = 100)
Heatmap(m_srt1.pct)

########################################### Avg Exp by cluster heatmap ########################################### 
m_srt1 <- m_srt
Idents(m_srt1) <- "celltype.group"

cluster.averages <- AverageExpression(m_srt, group.by = "celltype.group", features = gtc)
?AverageExpression
cluster.averages <- as.data.frame(cluster.averages)
colnames(cluster.averages) <- gsub(pattern = "RNA.", replacement = "", x = orig.levels)

cluster.averages$gene <- rownames(cluster.averages)
cluster.averages <- cluster.averages %>% gather(key='sample', value='value', -gene)
geom_tile()
ggplot(cluster.averages, aes(sample, gene)) + geom_tile(aes(fill=value)) + scale_fill_gradientn(colors = (c("#0200ad", "#fbfcbd", "#ff0000")))
#scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"), limits = c(0,2))

########################################### Avg Exp by cluster heatmap (seurat) ########################################### 
?AverageExpression
m_srt <- readRDS(file = "GFSPF_mito_ep_labelled_v2.rds")
Idents(m_srt) <- "group"
m_srt <- subset(m_srt, idents = "SPF")
Idents(m_srt) <- "celltype"
m_srt <- subset(m_srt, idents = c("ISC","TA1","TA2", "V1-V2","V3-V4","V5-V6"))
gtp <- c("Pycard","Casp1","Gsdmd","Il18","Il18r1")

cluster.averages1 <- AverageExpression(m_srt, group.by = "celltype", features = gtp, slot = "data",return.seurat = TRUE)
#cluster.averages1 <- AverageExpression(m_srt, add.ident ="celltype.group", features = gtc, return.seurat = TRUE)

Idents(cluster.averages1) <- factor(cluster.averages1$orig.ident, levels=c("ISC_GF","ISC_SPF", "TA1_GF","TA1_SPF", "TA2_GF", "TA2_SPF",
                                                                           "V1-V2_GF", "V1-V2_SPF","V3-V4_GF", "V3-V4_SPF","V5-V6_GF", "V5-V6_SPF", 
                                                                           "PC_GF", "PC_SPF","GC1_GF", "GC1_SPF","GC2_GF", "GC2_SPF", 
                                                                           "EE_GF", "EE_SPF","Tuft_GF","Tuft_SPF"))
Idents(cluster.averages1) <- factor(cluster.averages1$orig.ident, levels=c("RevSC_GF","RevSC_SPF","ISC_GF","ISC_SPF","EC1_GF","EC1_SPF",
                                    "EC2_GF","EC2_SPF","EC3_GF","EC3_SPF","GCPC_GF","GCPC_SPF","EE_GF","EE_SPF"))
Idents(cluster.averages1) <- factor(cluster.averages1$orig.ident, levels=c("ISC","TA1","TA2", "V1-V2","V3-V4","V5-V6","PC","GC1", "GC2","EE","Tuft"))

Idents(cluster.averages1) <- factor(cluster.averages1$orig.ident, levels=c("ISC_SPF","TA1_SPF","TA2_SPF", "V1-V2_SPF",
                                                                           "V3-V4_SPF","V5-V6_SPF","PC_SPF",
                                                                           "GC1_SPF", "GC2_SPF","EE_SPF","Tuft_SPF"))
Idents(cluster.averages1) <- factor(cluster.averages1$orig.ident, levels=c("ISC_SPF","TA1_SPF","TA2_SPF", "V1-V2_SPF",
                                                                           "V3-V4_SPF","V5-V6_SPF"))
Idents(cluster.averages1) <- factor(cluster.averages1$orig.ident, levels=c("RevSC","ISC","EC1", "EC2","EC3","GCPC","EE"))
cluster.averages1$orig.ident
?DoHeatmap
Idents(cluster.averages1)
png("SPF_mito_ep_Nate_HM_ECs.png", width = 500, height = 315, units = "px", res=100)
DoHeatmap(cluster.averages1, features = gtp, size = 5, draw.lines = FALSE, slot= "data",label = TRUE, group.by = "orig.ident",) +
  scale_fill_viridis()+ 
  theme(text = element_text(size = 16))
dev.off()
?DoHeatmap
#----------------seurat plot of individual cells------------------------------
Idents(m_srt1) <- "celltype.group"
DefaultAssay(m_srt) <- "SCT"
png("GFSPF_mito_ep_ISCTAPC_PRRs_HM_Allcells.png", width = 1000, height = 500, units = "px", res=100)
DoHeatmap(m_srt1, features = gtc, size=5,  cells = WhichCells(m_srt1, idents = c("RevSC_WT", "RevSC_Nod2KO"))) + 
  theme(text = element_text(size = 18)) +scale_fill_viridis()
dev.off()

#scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(m_srt, features = c("Ephb2"),
            order = FALSE, split.by = "group",
            cols = c("grey", "red"), pt.size = 2) & scale_colour_gradientn(colours = colours.for.plot)& theme(legend.position = "right")
dev.off
########################################### GSEA ########################################### 
Idents(m_srt) <- "celltype"
m_srt1 <- subset(m_srt, idents = c("GC1", "GC2"))
Idents(m_srt1) <- "group"
m_srt1.rank <- FindMarkers(m_srt1, ident.1 = "SPF", ident.2 = "GF", only.pos = FALSE,
                           logfc.threshold=0, min.pct=0, min.cells.group = 1, min.cells.feature = 1)
?FindMarkers
SPF_WT <- m_srt1.rank

m_srt1.rank <- SPF_WT
#m_srt1.rank$cluster <- ifelse(m_srt1.rank$avg_log2FC>0, "SPF", "GF")
m_srt1.rank$gene <- rownames(m_srt1.rank)

m_srt1.rank <- dplyr::filter(m_srt1.rank, avg_log2FC != 0)
count(m_srt1.rank)
m_srt1.rank <- dplyr::filter(m_srt1.rank, pct.1 & pct.2 != 0)
count(m_srt1.rank)

m_srt1.rank<- m_srt1.rank[, colnames(m_srt1.rank) != "pct.1"]
m_srt1.rank<- m_srt1.rank[, colnames(m_srt1.rank) != "pct.2"]
m_srt1.rank<- m_srt1.rank[, colnames(m_srt1.rank) != "p_val"]
m_srt1.rank<- m_srt1.rank[, colnames(m_srt1.rank) != "p_val_adj"]
#m_srt1.rank<- m_srt1.rank[, colnames(m_srt1.rank) != "cluster"]

geneList1 = m_srt1.rank[,1]
names(geneList1) = as.character(m_srt1.rank[,2])
geneList1 <- sort(geneList1, decreasing = TRUE)
#head(ranks)
#head(m_srt1.rank)
#head(geneList1)
ego3 <- gseGO(geneList     = geneList1,
              OrgDb        = org.Mm.eg.db,
              ont          = "BP",
              keyType      = "SYMBOL",
              scoreType = "std",
              eps          = 0,
              minGSSize    = 15,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

ridgeplot(ego3,  showCategory = 20,orderBy = "NES", label_format = 120)
write.csv(ego3, file = "GFSPF_mito_GC_CP_Gsea_GoBP.csv")

ridgeplot(ego3)
?ridgeplot
goplot(ego3)
dotplot(ego3)
cnetplot(ego3)

gseaplot(ego3, geneSetID = 2, by = "runningScore", title = ego3$Description[2])
ego3$d

goplot(kk2)

dotplot(kk2)
cnetplot(ego3)


#------------------GO gsea by entrez ID
m_srt1.EzID = bitr(m_srt1.rank$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
rownames(m_srt1.EzID) <- m_srt1.EzID$SYMBOL
m_srt1.EzID$avg_log2FC = m_srt1.rank$avg_log2FC[match(m_srt1.EzID$SYMBOL,m_srt1.rank$gene)]

geneList = m_srt1.EzID[,3]
names(geneList) = m_srt1.EzID[,2]
geneList <- sort(geneList, decreasing = TRUE)

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'mmu',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
ridgeplot(kk2)

wp2 <- gseWP(geneList, organism = "Mus musculus")
ridgeplot(wp2)

ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Mm.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

########################
ggplot(de, aes(avg_log2FC, -log10(p_val_adj))) +
  geom_point()  +
  geom_label_repel(data = de[ de$p_val_adj < 0.05 & de$avg_log2FC < 0.75, ], 
                   aes(label = gtc),
                   xlim = c(NA, -1), # <--- here
                   seed = 1) +
  geom_label_repel(data = de[ de$p_val_adj < 0.05 & de$avg_log2FC > 0.75, ], 
                   aes(label = gtc),
                   xlim = c(1, NA),  # <--- here
                   seed = 1) 
