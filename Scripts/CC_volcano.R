BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("AnnotationHub")
BiocManager::install("EnhancedVolcano")
BiocManager::install("ComplexHeatmap")
BiocManager::install("ggrepel")
BiocManager::install("enrichplot")

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
library(ggrepel)
library(enrichplot)

library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(digest)
library(cluster)
#GFSPF_mito_ep_labelled_v2
#EC_IR_GF_SPF_v2_label.rds
#GFSPF_IR_EC_SCT.rds
m_srt <- readRDS(file = "GFSPF_IR_EC_SCT.rds")
DimPlot(m_srt, reduction = "umap", group.by = "celltype",label = TRUE, label.size = 6)
m_srt1 <-m_srt
Idents(m_srt1) <- "celltype"
m_srt1 <- subset(m_srt, idents = c("RevSC"))
DimPlot(m_srt1, reduction = "umap", group.by = "celltype",label = TRUE, label.size = 6)
######################################################################
#volcano plots
######################################################################

Idents(m_srt1) <- "group"
Idents(m_srt1) <- "celltype.group"
m_srt1 <- PrepSCTFindMarkers(m_srt1)

#m_srt1.markers <- FindMarkers(m_srt1, ident.1 = "SPF", ident.2 = "GF", only.pos = FALSE,logfc.threshold=0.1, min.pct=0.1, min.cells.group = 1, min.cells.feature = 1)
#write.csv(m_srt1.markers, file = "GFSPF_mito_DE_Tuft.csv")

m_srt1.markers <- FindMarkers(m_srt1, ident.1 = "RevSC_SPF", ident.2 = "RevSC_GF", only.pos = FALSE, assay = "SCT",
                              logfc.threshold=0, min.pct=0, min.cells.group = 1, min.cells.feature = 1)

m_srt1.rank <- m_srt1.markers
m_srt1.markers$cluster <- ifelse(m_srt1.markers$avg_log2FC>0, "SPF", "GF")
m_srt1.markers$gene <- rownames(m_srt1.markers)
write.csv(m_srt1.markers, file = "GFSPF_IR_RevSC_DEallgenes.csv")

de <- as.data.frame(m_srt1.markers[,c("gene","avg_log2FC","p_val_adj", "cluster")]) 
#labels based on sig from fold exp and pvalue
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$avg_log2FC > 0.75 & de$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$avg_log2FC < -0.75 & de$p_val_adj < 0.05] <- "DOWN"

de$delabel[de$diffexpressed != "NO"] <- de$gene[de$diffexpressed != "NO"]
write.csv(de, file = "GFSPF_IR_RevSC_VolcanoLabels.csv")

gtc <- de$delabel
#gtc <- c("Mtor","Rptor","Mlst8")
png("GFSPF_IR_RevSC_Volcano1.png", width = 700, height = 800, units = "px", res=100)
EnhancedVolcano(de,
                lab = de$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = c(gtc),
                axisLabSize = 25,
                FCcutoff = 0.75,
                pCutoff = 0.05,
                labSize = 8,
                pointSize = 4,
                col = c("grey30", "grey30", "grey30", "red2"),
                drawConnectors = TRUE,
                maxoverlapsConnectors = 10)
dev.off()
#drawConnectors = TRUE
#directionConnectors = "y"
#maxoverlapsConnectors = Inf
?EnhancedVolcano
######################################################################
#GSEA
######################################################################

#Idents(m_srt) <- "celltype"
#m_srt1 <- subset(m_srt, idents = c("EC1", "EC2"))
#Idents(m_srt1) <- "group"
#m_srt1.rank <- FindMarkers(m_srt1, ident.1 = "SPF", ident.2 = "GF", only.pos = FALSE, logfc.threshold=0, min.pct=0, min.cells.group = 1, min.cells.feature = 1)

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
ego1 <- gseGO(geneList     = geneList1,
              OrgDb        = org.Mm.eg.db,
              ont          = "BP",
              keyType      = "SYMBOL",
              scoreType = "std",
              eps          = 0,
              minGSSize    = 15,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

ego1 <- pairwise_termsim(ego1)
ego1 <- simplify(ego1, cutoff=0.7, by="p.adjust", select_fun=min)

png("GFSPF_IR_RevSC_GoBP_wnt.png", width = 800, height = 500, units = "px", res=100)
ridgeplot(ego1,  showCategory = 12, orderBy = "NES", label_format = 45)+
  theme(axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.key.size = unit(5, 'mm'))
dev.off()
#?ridgeplot

write.csv(ego1, file = "GFSPF_IR_RevSC_Gsea_GoBP.csv")
gseaplot(ego1, geneSetID = 15, by = "runningScore", title = ego1$Description[1])

ego2 <- gseGO(geneList     = geneList1,
              OrgDb        = org.Mm.eg.db,
              ont          = "CC",
              keyType      = "SYMBOL",
              scoreType = "std",
              eps          = 0,
              minGSSize    = 15,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

png("GFSPF_IR_RevSC_Gsea_GoCC.png", width = 800, height = 400, units = "px", res=100)
ridgeplot(ego2,  showCategory = 10,orderBy = "NES", label_format = 60) 
dev.off()
write.csv(ego2, file = "GFSPF_IR_RevSC_Gsea_GoCC.csv")

######################################################################
#Compare Cluster
######################################################################

m_srt1.CC <- FindMarkers(m_srt1, ident.1 = "SPF", ident.2 = "GF", only.pos = FALSE, verbose = FALSE)
m_srt1.CC  %>% top_n(n = 2, wt = avg_log2FC)
m_srt1.CC$cluster <- ifelse(m_srt1.markers$avg_log2FC>0, "SPF", "GF")

pval <- subset(m_srt1.CC, rowSums(m_srt1.markers[5] < 0.05) > 0)
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

png("GFSPF_IR_RevSC_CompCluster_eGO_dotplot.png", width = 600, height = 600, units = "px", res=100)
dotplot(GOclusterplot, label_format = 40)+theme(axis.text.y = element_text(size = 15))
dev.off()

write.csv(GOclusterplot, file = "GFSPF_IR_RevSC_CompCluster_eGO.csv")

#enrichWPclusterplot <- compareCluster(geneCluster = genelist1, fun = "enrichWP", organism = "Mus musculus")
#write.csv(enrichWPclusterplot, file = "GFSPF_mito_AllCells_WPcluster.csv")
#dotplot(enrichWPclusterplot)

#KEGGclusterplot <- compareCluster(geneCluster = genelist1, fun = "enrichKEGG")
#dotplot(KEGGclusterplot)