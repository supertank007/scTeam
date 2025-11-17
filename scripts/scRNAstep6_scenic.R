##可视化
rm(list=ls())
setwd("/DATA/")
source("scoption.cfg")
library(Seurat)
library(BiocManager)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)
load("STEP3.RData")
#### 1.提取 out_SCENIC.loom 信息
loom <- open_loom('./scenic/sample_SCENIC.loom') 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])
embeddings <- get_embeddings(loom)  
close_loom(loom)
rownames(regulonAUC)
names(regulons)
sub_regulonAUC <- regulonAUC[,match(colnames(Cell_Clusters),colnames(regulonAUC))]
dim(sub_regulonAUC)

#确认是否一致
identical(colnames(sub_regulonAUC), colnames(Cell_Clusters))
cellClusters <- data.frame(row.names = colnames(Cell_Clusters), 
                           seurat_clusters = as.character(Cell_Clusters$seurat_clusters))
cellTypes <- data.frame(row.names = colnames(Cell_Clusters), 
                        celltype = Cell_Clusters$celltype)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 

#保存一下
save(sub_regulonAUC,cellTypes,cellClusters,seurat.data,
     file = 'for_rss_and_visual.Rdata')


### 4.1. TF活性均值
# 看看不同单细胞亚群的转录因子活性平均值
# Split the cells by cluster:
selectedResolution <- "celltype" # select resolution
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])

# 去除extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
pdf(file.path("scenic/", "Heatmap.pdf"),width = 9,height = 25)
Heatmap(
  regulonActivity_byGroup_Scaled,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
dev.off()

.H <- function(pVect){
  pVect <- pVect[pVect>0] # /sum(pVect) ??
  - sum(pVect * log2(pVect))
}

# Jensen-Shannon Divergence (JSD)
calcJSD <- function(pRegulon, pCellType)
{
  (.H((pRegulon+pCellType)/2)) - ((.H(pRegulon)+.H(pCellType))/2)
}

# Regulon specificity score (RSS)
.calcRSS.oneRegulon <- function(pRegulon, pCellType)
{
  jsd <- calcJSD(pRegulon, pCellType)
  1 - sqrt(jsd)
}

calcRSS <- function(AUC, cellAnnotation, cellTypes=NULL)
{
  if(any(is.na(cellAnnotation))) stop("NAs in annotation")
  if(any(class(AUC)=="aucellResults")) AUC <- getAUC(AUC)
  normAUC <- AUC/rowSums(AUC)
  if(is.null(cellTypes)) cellTypes <- unique(cellAnnotation)
  # 
  ctapply <- lapply
  if(require('BiocParallel')) ctapply <- bplapply
  
  rss <- ctapply(cellTypes, function(thisType)
    sapply(rownames(normAUC), function(thisRegulon)
    {
      pRegulon <- normAUC[thisRegulon,]
      pCellType <- as.numeric(cellAnnotation==thisType)
      pCellType <- pCellType/sum(pCellType)
      .calcRSS.oneRegulon(pRegulon, pCellType)
    })
  )
  rss <- do.call(cbind, rss)
  colnames(rss) <- cellTypes
  return(rss)
}

rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellTypes[colnames(sub_regulonAUC), selectedResolution]) 
rss=na.omit(rss) 

plotRSS <- function(rss, labelsToDiscard=NULL, zThreshold=1,
                    cluster_columns=FALSE, order_rows=TRUE, thr=0.01, varName="cellType",
                    col.low="grey90", col.mid="darkolivegreen3", col.high="darkgreen",
                    revCol=FALSE, verbose=TRUE)
{
  varSize="RSS"
  varCol="Z"
  if(revCol) {
    varSize="Z"
    varCol="RSS"
  }

    rssNorm <- scale(rss) # scale the full matrix...
    rssNorm <- rssNorm[,which(!colnames(rssNorm) %in% labelsToDiscard)] # remove after calculating...
    rssNorm[rssNorm < 0] <- 0
    
    ## to get row order (easier...)
    rssSubset <- rssNorm
    if(!is.null(zThreshold)) rssSubset[rssSubset < zThreshold] <- 0
    tmp <- .plotRSS_heatmap(rssSubset, thr=thr, cluster_columns=cluster_columns, order_rows=order_rows, verbose=verbose)
    rowOrder <- rev(tmp@row_names_param$labels)
    rm(tmp)  
    ## Dotplot
    rss.df <- reshape2::melt(rss)
    head(rss.df)
    colnames(rss.df) <- c("Topic", varName, "RSS")
    rssNorm.df <- reshape2::melt(rssNorm)
    colnames(rssNorm.df) <- c("Topic", varName, "Z")
    rss.df <- base::merge(rss.df, rssNorm.df)
    
    rss.df <- rss.df[which(!rss.df[,varName] %in% labelsToDiscard),] # remove after calculating...
    if(nrow(rss.df)<2) stop("Insufficient rows left to plot RSS.")
    
    rss.df <- rss.df[which(rss.df$Topic %in% rowOrder),]
    rss.df[,"Topic"] <- factor(rss.df[,"Topic"], levels=rowOrder)
    p <- dotHeatmap(rss.df, 
                    var.x=varName, var.y="Topic", 
                    var.size=varSize, min.size=.5, max.size=5,
                    var.col=varCol, col.low=col.low, col.mid=col.mid, col.high=col.high)
    
    invisible(list(plot=p, df=rss.df, rowOrder=rowOrder))
}


.plotRSS_heatmap <- plotRSS_heatmap <- function(rss, thr=NULL, row_names_gp=gpar(fontsize=5), order_rows=TRUE, cluster_rows=FALSE, name="RSS", verbose=TRUE, ...)
{
  if(is.null(thr)) thr <- signif(quantile(rss, p=.97),2)
  
  library(ComplexHeatmap)
  rssSubset <- rss[rowSums(rss > thr)>0,]
  rssSubset <- rssSubset[,colSums(rssSubset > thr)>0]
  
  if(verbose) message("Showing regulons and cell types with any RSS > ", thr, " (dim: ", nrow(rssSubset), "x", ncol(rssSubset),")")
  
  if(order_rows)
  {
    maxVal <- apply(rssSubset, 1, which.max)
    rss_ordered <- rssSubset[0,]
    for(i in 1:ncol(rssSubset))
    {
      tmp <- rssSubset[which(maxVal==i),,drop=F]
      tmp <- tmp[order(tmp[,i], decreasing=FALSE),,drop=F]
      rss_ordered <- rbind(rss_ordered, tmp)
    }
    rssSubset <- rss_ordered
    cluster_rows=FALSE
  }
  
  Heatmap(rssSubset, name=name, row_names_gp=row_names_gp, cluster_rows=cluster_rows, ...)
} 

dotHeatmap <- function (enrichmentDf,
                        var.x="Topic", var.y="ID", 
                        var.col="FC", col.low="dodgerblue", col.mid="floralwhite", col.high="brown1", 
                        var.size="p.adjust", min.size=1, max.size=8,
                        ...)
{
  require(data.table)
  require(ggplot2)
  
  colorPal <- grDevices::colorRampPalette(c(col.low, col.mid, col.high))
  p <- ggplot(data=enrichmentDf, mapping=aes_string(x=var.x, y=var.y)) + 
    geom_point(mapping=aes_string(size=var.size, color=var.col)) +
    scale_radius(range=c(min.size, max.size)) +
    scale_colour_gradientn(colors=colorPal(10)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y=element_blank(), 
          axis.text.x=element_text(angle=90, hjust=1),
          ...)
  return(p)
}

rssPlot <- plotRSS(rss)
pdf(file.path("scenic/", "rss.pdf"),width = 9,height = 25)
rssPlot$plot
dev.off()

rss=regulonActivity_byGroup_Scaled
head(rss)
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss),
                 cluster =   colnames(rss)[i],
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster) 
n = rss[top5$path,] 
#rownames(rowcn) = rownames(n)
pdf(file.path("scenic/", "top5.pdf"),width = 9,height = 15)
pheatmap(n,
         annotation_row = rowcn,
         show_rownames = T)
dev.off()
