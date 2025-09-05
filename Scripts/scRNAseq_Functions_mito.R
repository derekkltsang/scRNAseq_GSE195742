
#Seurat dependency
# ---------------- Differential Gene Exp -------------------------------
diff_gene_exp <- function(object, cluster_ID, ident_1, ident_2, name){
DefaultAssay(object) <- "RNA"
Idents(object) <- "celltype"

g <- subset(object, idents = cluster_ID)
Idents(g) <- "group"

avg.g <- log1p(AverageExpression(g, verbose = FALSE)$RNA)
avg.g$gene <- rownames(avg.g)

Idents(object) <- "celltype.group"
response <- FindMarkers(object, ident.1 = ident_1, ident.2 = ident_2, verbose = FALSE)
head(response, n = 15)

csvname <- paste(c(name, ".csv"), collapse = "")
write.csv(response, csvname)
gtp <- rownames(response[(response$p_val_adj < 0.001),0])

pngname <- paste(c(name, ".png"), collapse = "")

png(pngname, width = 500, height = 500, units = "px")
p <- ggplot(avg.g, aes(GF, SPF)) + geom_point() + ggtitle(cluster_ID)
p <- LabelPoints(plot = p, points = gtp)
print(p)
dev.off()
}


diff_gene_exp <- function(object, cluster_ID, ident_1, ident_2, name){
  DefaultAssay(object) <- "RNA"
  Idents(object) <- "celltype"
  
  g <- subset(object, idents = cluster_ID)
  Idents(g) <- "group"
  
  avg.g <- as.data.frame(log1p(AverageExpression(g, group.by= "group", verbose = FALSE)$RNA))
  #avg.g$gene <- rownames(avg.g)
  
  Idents(object) <- "group"
  response <- FindMarkers(object, ident.1 = ident_1, ident.2 = ident_2, verbose = FALSE)
  head(response, n = 15)
  
  csvname <- paste(c(name, cluster_ID, ".csv"), collapse = "")
  write.csv(response, csvname)
  gtp <- rownames(response[(response$p_val_adj < 0.001),0])
  
  pngname <- paste(c(name, cluster_ID, ".png"), collapse = "")
  
  png(pngname, width = 500, height = 500, units = "px")
  p <- ggplot(avg.g, aes(GF, SPF)) + geom_point() + ggtitle(cluster_ID)
  p <- LabelPoints(plot = p, points = gtp)
  print(p)
  dev.off()
}
# -------------------------------------------------------------------



# ---------------- Differential Gene Exp -------------------------------
topGO_diff <- function (cluster_ID, ident_1, ident_2, csv_name) {
  
Idents(m_srt) <- "celltype"
cluster <- subset(m_srt, idents = cluster_ID)

Idents(cluster) <- "group"
DE <- FindMarkers(cluster, ident.1 = ident_1, ident.2 = ident_2, verbose = FALSE, only.pos = TRUE)

?FindMarkers

expr <- cluster@assays$RNA@data
n.gt.0 <- apply(expr, 1, function(x)length(which(x > 0)))
all.genes <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.10)] #find genes expressed in more than 10% of cells

ag <- DE$p_val
ag <- rownames(DE)
  
ag2 <- ifelse(all.genes %in% ag, 1, 0)
names(ag2) <- all.genes
  
  
GOdata <- new("topGOdata",
                ontology = "BP", # use biological process ontology
                allGenes = ag2,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.org, 
                mapping = "org.Mm.eg.db", 
                ID = "symbol")
  
# Test for enrichment using Fisher's Exact Test
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

topgotable <- GenTable(GOdata, Fisher = resultFisher, topNodes = 50, numChar = 60)

write.csv(topgotable, csv_name)
}
# -------------------------------------------------------------------

colours.for.plot <- c("#d3d3d3", "#ffeda0" ,"#fed976","#feb24c",
                      "#fd8d3c", "#fc4e2a", "#e31a1c" ,"#bd0026","#800026")



FPlot_split <- function (genetoplot, tf_order, name, maxcutoff) {
  
  gtp1 <- FeaturePlot(m_srt1, genetoplot,
                      min.cutoff = 0, max.cutoff = maxcutoff, split.by = "group",
                      order = tf_order, pt.size = 2.5) + scale_colour_gradientn(colours = colours.for.plot) + theme(legend.position = "right")
  
  gtp2 <- FeaturePlot(m_srt2, features = genetoplot,
                      min.cutoff = 0, max.cutoff = maxcutoff, split.by = "group",
                      order = tf_order, pt.size = 2.5) + scale_colour_gradientn(colours = colours.for.plot) + theme(legend.position = "right")
  
  pngname <- paste(c(name, genetoplot, ".png"), collapse = "")
  
  png(pngname, width = 1000, height = 500, units = "px", pointsize= 20, res = 100)
  
  print(plot_grid(gtp1, gtp2))
  
  dev.off()
  
}

# -------------------------------------------------------------------
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    (sum(counts[genes,]>0)/ncells)*100
  }else{return(NA)}
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

#raw counts (not normalized by cell counts)
calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)
  }else{return(NA)}
}
# -------------------------------------------------------------------




`plot_seperate_features <- function(object, gene, limit = NULL,ident=NULL,pt.size=0.5,revers.order = F,numcol=NULL, colors =c('#DCD9DE','#870052'), ...){
  
  if(is.null(ident)){
    Idents(object) <- 'combined'
  }
  else{
    Idents(object) <- object[[ident]]
  }
  p <- vector("list", length = (length(gene)*length(levels(Idents(object)))))
  for(g in 1:(length(gene))){
    if(is.null(limit)==T){
      for(i in 1:length(levels(Idents(object)))){
        j <- i+((g-1)*length(levels(Idents(object))))
        p[[j]] <- FeaturePlot(object, cells = WhichCells(object,idents = levels(Idents(object))[i]), features = gene[g],order = T,cols = colors, pt.size = pt.size, sort.cell = T,...) +
          scale_color_gradient(low = colors[1], high = colors[2], limits = c(0,ceiling(max(object@assays$RNA@data[gene[g],]))),oob = squish) +
          ggtitle(paste(levels(Idents(object))[i],gene[g],sep=", "))
      }
    }
    else{
      for(i in 1:length(levels(Idents(object)))){
        j <- i+((g-1)*length(levels(Idents(object))))
        p[[j]] <- FeaturePlot(object, cells = WhichCells(object,idents = levels(Idents(object))[i]), features = gene[g],order = T,cols = colors, pt.size = pt.size, sort.cell = T,...) +
          scale_color_gradient(low = colors[1], high = colors[2], limits = limit,oob = squish) +
          ggtitle(paste(levels(Idents(object))[i],gene[g],sep=", "))
      }
    }
  }
  #correct axis in all plots
  y_axis_max <- 1:length(p)
  y_axis_min <- 1:length(p)
  x_axis_max <- 1:length(p)
  x_axis_min <- 1:length(p)
  for(i in 1:length(p)){
    x_axis_max[i] <- max(layer_scales(p[[i]])$x$range$range)
    x_axis_min[i] <- min(layer_scales(p[[i]])$x$range$range)
    y_axis_max[i] <- max(layer_scales(p[[i]])$y$range$range)
    y_axis_min[i] <- min(layer_scales(p[[i]])$y$range$range)
  }
  for(i in 1:length(p)){
    p[[i]] <- p[[i]]+coord_equal(xlim=c(min(x_axis_min),max(x_axis_max)),ylim=c(min(y_axis_min),max(y_axis_max)))
  }
  if(revers.order==T) p <- rev(p)
  gridExtra::grid.arrange(grobs = p,ncol=numcol)
}`