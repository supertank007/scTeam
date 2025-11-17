#环境清理及载入包
rm(list = ls())
library(Seurat)
library(cowplot)
library(dplyr)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(DoubletFinder)
library(clustree)


setwd(Sys.getenv("PWD"))
source("scoption.cfg")
#读入数据单细胞数据
dir_list <- list.dirs()
dir <- sub("^\\./", "", dir_list)
dir <- dir[dir != "."]
dir <- dir[dir != "results"]
dir <- dir[dir != "results/Step1results"]
dir <- dir[dir != "results/Step2results"]
dir <- dir[dir != "results/markers"]
dir <- dir[dir != "results/SingleR"]
dir <- dir[dir != "scenic"]
dir <- dir[dir != "out"]

#读入
counts <- Read10X(data.dir = dir)#根据文件不同修改读入方式

#创建项目
DATA.list <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  DATA.list[[i]] <- CreateSeuratObject(counts, min.cells = min_cells, min.features = min_features)
}

#质控
dir.create("results")
dir.create("results/Step1results")
if (genome == "mouse") {
  for(i in 1:length(dir)){
  DATA.list[[i]][["percent.mt"]] <- PercentageFeatureSet(DATA.list[[i]], pattern = "^mt-")
  VlnPlot(DATA.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  file_name <- paste0("VlnPlot_", dir[i], ".pdf")
  ggsave(file.path("results/Step1results", file_name), device = "pdf")
  DATA.list[[i]] <- subset(DATA.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < max_features & percent.mt < percent_mt)
  }
} else if (genome == "human") {
  for(i in 1:length(dir)){
  DATA.list[[i]][["percent.mt"]] <- PercentageFeatureSet(DATA.list[[i]], pattern = "^MT-")
  VlnPlot(DATA.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  file_name <- paste0("VlnPlot_", dir[i], ".pdf")
  ggsave(file.path("results/Step1results", file_name), device = "pdf")
  DATA.list[[i]] <- subset(DATA.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < max_features & percent.mt < percent_mt)
  }
} else {
  system("echo 'Genome information not set'", intern = TRUE)
  stop("Please verify and set the genome. Exiting R.")
}



#赋予组别
for(i in 1:length(dir)){DATA.list[[i]]$id <- dir[i]}

#联合分析
#标准化寻找高变基因
DATA.list <- lapply(X = DATA.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#去批次整合
DATA.anchors <- FindIntegrationAnchors(object.list = DATA.list, dims = 1:20)
DATA.combined <- IntegrateData(anchorset = DATA.anchors, dims = 1:20)

#进入分析-选择分群分辨率
DefaultAssay(DATA.combined) <- "integrated"
DATA.combined <- ScaleData(DATA.combined, verbose = FALSE)
DATA.combined <- RunPCA(DATA.combined, features = VariableFeatures(object = DATA.combined))
DATA.combined <- RunUMAP(DATA.combined, dims = 1:30)
DATA.combined <- RunTSNE(DATA.combined, dims = 1:30)
DATA.combined <- FindNeighbors(DATA.combined, reduction = "pca", dims = 1:30)
#绘制分群脉络图
DATA.combined <- FindClusters(DATA.combined, resolution = c(seq(.2,2,.3)))
clust <- clustree(DATA.combined@meta.data, prefix = "integrated_snn_res.")
ggsave(file.path("results/Step1results", "clusttree.pdf"), plot = clust, width = 8, height = 12)
PC1_2 <- VizDimLoadings(DATA.combined, dims = 1:2, reduction = "pca")
ggsave(file.path("results/Step1results", "PC.pdf"), plot = PC1_2, width = 6, height = 6)
PCA <- DimPlot(DATA.combined, reduction = "pca")
ggsave(file.path("results/Step1results", "PCA.pdf"), plot = PCA, width = 9, height = 6)
pdf(file.path("results/Step1results", "heatmap_plot.pdf"))
DimHeatmap(DATA.combined, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()
DATA.combined <- JackStraw(DATA.combined, num.replicate = 100)
DATA.combined <- ScoreJackStraw(DATA.combined, dims = 1:20)
JackStraw <- JackStrawPlot(DATA.combined, dims = 1:15)
ggsave(file.path("results/Step1results", "JackStraw.pdf"), plot = JackStraw, width = 9, height = 18)
Elbow <- ElbowPlot(DATA.combined)
ggsave(file.path("results/Step1results", "Elbow.pdf"), plot = Elbow, width = 9, height = 6)
rm(Elbow)
rm(JackStraw)
rm(PCA)
rm(PC1_2)
rm(clust)
rm(DATA.anchors)
rm(counts)
rm(DATA.list)
save.image(file = "STEP1.RData")
