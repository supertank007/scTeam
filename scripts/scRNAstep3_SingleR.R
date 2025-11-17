#SingleR注释
setwd("/DATA")
source("scoption.cfg")
load("STEP2.RData")
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
#(mouseRNA,blueRNA,DatabaseImmuneCell，HumanPrimaryCellAtlas，ImmGenData ，MonacoImmune，BlueprintEncode，NovershternHematopoietic)


load("/DATABANK/ImmGenData.Rdata")
load("/DATABANK/DatabaseImmuneCell.Rdata")
load("/DATABANK/MouseRNAseqData.Rdata")
load("/DATABANK/BlueprintEncode.Rdata")
load("/DATABANK/HumanPrimaryCellAtlas.Rdata")
load("/DATABANK/ImmGenData.Rdata")
load("/DATABANK/NovershternHematopoietic.Rdata")
load("/DATABANK/MonacoImmune.Rdata")
dir.create("results/SingleR")

if (genome == "mouse") {
  sce_for_SingleR <- GetAssayData(Cell_Clusters, slot="data")
  pred.mouseRNA<- SingleR(test = sce_for_SingleR, 
                        ref = mouseRNA, 
                        clusters = Cell_Clusters$seurat_clusters,
                        labels = mouseRNA$label.main)
  table(pred.mouseRNA$labels)
  Cell_Clusters$singleR_cluster <- pred.mouseRNA$labels[match(Cell_Clusters$seurat_clusters, rownames(pred.mouseRNA))]
  table(Cell_Clusters$singleR_cluster)
  SingleRanntsne <- DimPlot(Cell_Clusters, group.by = c("singleR_cluster"),reduction = "tsne")
  SingleRannumap <- DimPlot(Cell_Clusters, group.by = c("singleR_cluster"),reduction = "umap")
  SingleRanntsne_split <- DimPlot(Cell_Clusters, group.by = c("singleR_cluster"),reduction = "tsne", split.by  = "Group")
  SingleRannumap_split <- DimPlot(Cell_Clusters, group.by = c("singleR_cluster"),reduction = "umap", split.by  = "Group")
  ggsave(file.path("results/SingleR","SingleR_tsne.pdf"),plot = SingleRanntsne)
  ggsave(file.path("results/SingleR","SingleR_umap.pdf"),plot = SingleRannumap)
  ggsave(file.path("results/SingleR","SingleR_tsne_split.pdf"),plot = SingleRanntsne_split)
  ggsave(file.path("results/SingleR","SingleR_umap_split.pdf"),plot = SingleRannumap_split)
} else if (genome == "human") {
  sce_for_SingleR <- GetAssayData(Cell_Clusters, slot="data")
  pred.humanRNA<- SingleR(test = sce_for_SingleR, 
                        ref = HumanPrimaryCellAtlas, 
                        clusters = Cell_Clusters$seurat_clusters,
                        labels = HumanPrimaryCellAtlas$label.main)
  table(pred.humanRNA$labels)
  Cell_Clusters$singleR_cluster <- pred.humanRNA$labels[match(Cell_Clusters$seurat_clusters, rownames(pred.humanRNA))]
  table(Cell_Clusters$singleR_cluster)
  SingleRanntsne <- DimPlot(Cell_Clusters, group.by = c("singleR_cluster"),reduction = "tsne")
  SingleRannumap <- DimPlot(Cell_Clusters, group.by = c("singleR_cluster"),reduction = "umap")
  SingleRanntsne_split <- DimPlot(Cell_Clusters, group.by = c("singleR_cluster"),reduction = "tsne", split.by  = "Group")
  SingleRannumap_split <- DimPlot(Cell_Clusters, group.by = c("singleR_cluster"),reduction = "umap", split.by  = "Group")
  ggsave(file.path("results/SingleR","SingleR_tsne.pdf"),plot = SingleRanntsne)
  ggsave(file.path("results/SingleR","SingleR_umap.pdf"),plot = SingleRannumap)
  ggsave(file.path("results/SingleR","SingleR_tsne_split.pdf"),plot = SingleRanntsne_split)
  ggsave(file.path("results/SingleR","SingleR_umap_split.pdf"),plot = SingleRannumap_split)
} else {
  system("echo 'Genome information not set'", intern = TRUE)
  stop("Please verify and set the genome. Exiting R.")
}

####自动注释完成
#重新注释
