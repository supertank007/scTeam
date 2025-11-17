setwd(Sys.getenv("PWD"))
source("scoption.cfg")

library(Seurat)
library(cowplot)
library(dplyr)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(clustree)
#确定好参数后设置好分辨率进行分群
#设置分组
if (automatic_resolution == 1) {
  load("Resolution.RData")
} else if (manual_resolution == 0) {
  load("STEP1.RData")
  Cell_Clusters <- FindClusters(DATA.combined, resolution = resolution)
  rm(DATA.combined)
}

Cell_Clusters$Group <- ""
Cell_Clusters$Group[Cell_Clusters$id %in% Group1] <- Group1name
Cell_Clusters$Group[Cell_Clusters$id %in% Group2] <- Group2name
Cell_Clusters$Group[Cell_Clusters$id %in% Group3] <- Group3name
Cell_Clusters$Group[Cell_Clusters$id %in% Group4] <- Group4name
#降维分析
plot1 <- DimPlot(Cell_Clusters, reduction = "umap", group.by = "id")
plot2 <- DimPlot(Cell_Clusters, reduction = "pca", label = TRUE)
plot3 <- DimPlot(Cell_Clusters, reduction = "tsne", label = TRUE)
plot4 <- DimPlot(Cell_Clusters, reduction = "umap", label = TRUE)
dir.create("results/Step2results")
ggsave(file.path("results/Step2results","Oringeplot.pdf"), plot = plot1+plot2+plot3+plot4, width = 25, height = 25)
ggsave(file.path("results/Step2results","umap_group.pdf"), plot = plot1, width = 9, height = 9)
ggsave(file.path("results/Step2results","pca_or.pdf"), plot = plot2, width = 9, height = 9)
ggsave(file.path("results/Step2results","tsne_or.pdf"), plot = plot3, width = 9, height = 9)
ggsave(file.path("results/Step2results","umap_or.pdf"), plot = plot4, width = 9, height = 9)
#SPLIT
plot1 <- DimPlot(Cell_Clusters, reduction = "umap", split.by  = "id")
plot2 <- DimPlot(Cell_Clusters, reduction = "pca", split.by  = "id" ,label = TRUE)
plot3 <- DimPlot(Cell_Clusters, reduction = "tsne", split.by  = "id" ,label = TRUE)
plot4 <- DimPlot(Cell_Clusters, reduction = "umap", split.by  = "id",label = TRUE)
ggsave(file.path("results/Step2results","splitplot.pdf"), plot = plot1+plot2+plot3+plot4, width = 25, height = 25)
ggsave(file.path("results/Step2results","umap_split.pdf"), plot = plot1, width = 9, height = 9)
ggsave(file.path("results/Step2results","pca_split.pdf"), plot = plot2, width = 9, height = 9)
ggsave(file.path("results/Step2results","tsne_split.pdf"), plot = plot3, width = 9, height = 9)
ggsave(file.path("results/Step2results","umap_split_lable.pdf"), plot = plot4, width = 9, height = 9)


plot1 <- DimPlot(Cell_Clusters, reduction = "umap", split.by  = "Group")
plot2 <- DimPlot(Cell_Clusters, reduction = "pca", split.by  = "Group" ,label = TRUE)
plot3 <- DimPlot(Cell_Clusters, reduction = "tsne", split.by  = "Group" ,label = TRUE)
plot4 <- DimPlot(Cell_Clusters, reduction = "umap", split.by  = "Group",label = TRUE)
ggsave(file.path("results/Step2results","splitplot_Group.pdf"), plot = plot1+plot2+plot3+plot4, width = 25, height = 25)
ggsave(file.path("results/Step2results","umap_split_Group.pdf"), plot = plot1, width = 9, height = 9)
ggsave(file.path("results/Step2results","pca_split_Group.pdf"), plot = plot2, width = 9, height = 9)
ggsave(file.path("results/Step2results","tsne_split_Group.pdf"), plot = plot3, width = 9, height = 9)
ggsave(file.path("results/Step2results","umap_split_lable_Group.pdf"), plot = plot4, width = 9, height = 9)

#计算细胞比例
table(Cell_Clusters$Group)#查看各组细胞数
prop.table(table(Idents(Cell_Clusters)))
table(Idents(Cell_Clusters), Cell_Clusters$Group)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(Cell_Clusters), Cell_Clusters$Group), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
cellrat <- ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave(file.path("results/Step2results","cellratio.pdf"), plot = cellrat, width = 8, height = 5)

#转录组分析
dir.create("results/markers")
DefaultAssay(Cell_Clusters) <- "RNA"
Cell_Clusters <- NormalizeData(Cell_Clusters, verbose=F)
Cell_Clusters.markers <- FindAllMarkers(Cell_Clusters, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
Cell_Clusters.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- Cell_Clusters.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- Cell_Clusters.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
Cell_Clusters@assays$RNA@scale.data <- scale(Cell_Clusters@assays$RNA@data, scale = TRUE)
pdf(file.path("results/markers", "top10heatmapgene.pdf"))
DoHeatmap(Cell_Clusters, features = top10$gene) + NoLegend()
dev.off()
write.csv(Cell_Clusters.markers, file.path("results/markers", "all.markers.csv"),row.names = F)
write.csv(top10, file.path("results/markers","top10_markers.csv"),row.names = F)
write.csv(top2, file.path("results/markers","top2_markers.csv"),row.names = F)
rm(top2)
rm(top10)
rm(cellrat)
rm(Cell_Clusters.markers)
rm(Cellratio)
rm(plot1)
rm(plot2)
rm(plot3)
rm(plot4)
rm(allcolour)
save.image(file = "STEP2.RData")