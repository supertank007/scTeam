#人工注释
setwd(Sys.getenv("PWD"))
load("STEP2.RData")
library(Seurat)
library(ggplot2)
source("scoption.cfg")
new.cluster.ids <- new_cluster_ids

names(new.cluster.ids) <- levels(Cell_Clusters)
Cell_Clusters <- RenameIdents(Cell_Clusters, new.cluster.ids)
Cell_Clusters$celltype <- new.cluster.ids[Idents(Cell_Clusters)]
annotsne <- DimPlot(Cell_Clusters, reduction = 'tsne', split.by = c("Group"), label = TRUE, pt.size = 0.5)
ggsave(file.path("results/","annotsne.pdf"),plot = annotsne,width = 10,height = 5)
annoumap <- DimPlot(Cell_Clusters, reduction = 'umap', split.by = c("Group"), label = TRUE, pt.size = 0.5)
ggsave(file.path("results/","annoumap.pdf"),plot = annoumap,width = 10,height = 5)

#anno <- DimPlot(Cell_Clusters, reduction = 'tsne', split.by = c("Group"),
        #label = TRUE, pt.size = 0.5, cols = c("#66C5CC","#F6CF71","#F89C74","#DCB0F2","#87C55F","#9EB9F3","#FE88B1","#C9DB74","#8BE0A4","#B6992D","#F1CE63","#499894")) + NoLegend()

#计算细胞比例
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
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 0.5, linetype="solid"))
ggsave(file.path("results/","cellratio_anno.pdf"), plot = cellrat, width = 8, height = 5)


#####marker展示

DefaultAssay(Cell_Clusters) <- "RNA"
markerfeaturetsne <- FeaturePlot(Cell_Clusters, features = markers, pt.size = 0, combine = T,reduction = "tsne")
markerfeaturetumap <- FeaturePlot(Cell_Clusters, features = markers, pt.size = 0, combine = T,reduction = "umap")
markervln <- VlnPlot(Cell_Clusters, features = markers, pt.size = 0, combine = T)
ggsave(file.path("results/","markerfeaturetsne.pdf"), plot = markerfeaturetsne, width = 16, height = 16)
ggsave(file.path("results/","markerfeaturetumap.pdf"), plot = markerfeaturetumap, width = 16, height = 16)
ggsave(file.path("results/","markervln.pdf"), plot = markervln, width = 16, height = 16)
write.csv(t(as.matrix(Cell_Clusters@assays$RNA@counts)),file = "scTF.csv")
rm(annotsne)
rm(annoumap)
rm(Cellratio)
rm(cellrat)
rm(markerfeaturetsne)
rm(markerfeaturetumap)
rm(markervln)
save.image(file = "STEP3.RData")
