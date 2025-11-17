library(tidyverse)
library(foreach)
library(doParallel)
library(dplyr)
library(Seurat)
library(cluster)
library(parallel)

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
Cell_cluster <- IntegrateData(anchorset = DATA.anchors, dims = 1:20)
DefaultAssay(Cell_cluster) <- "integrated"
Cell_cluster <- ScaleData(Cell_cluster)
Cell_cluster <- RunPCA(Cell_cluster,dims =1:50)

#一阶段分辨率选择
var_explained <- Cell_cluster[["pca"]]@stdev^2
var_explained_ratio <- var_explained / sum(var_explained)
cumulative_variance <- cumsum(var_explained_ratio)
threshold <- 0.9
optimal_dim <- which(cumulative_variance >= threshold)[1]
Cell_cluster <- RunPCA(Cell_cluster,dims = 1:optimal_dim)

#二阶段Kparam选择
find_optimal_k <- function(Cell_cluster, dims = 1:30, k_range = seq(10, 50, by = 5)) {
  results <- data.frame(k.param = integer(), mean_connectivity = numeric())
  
  for (k in k_range) {
    # 显式指定图名避免歧义
    graph_name <- paste0(DefaultAssay(Cell_cluster), "_snn")
    
    Cell_cluster <- FindNeighbors(Cell_cluster, k.param = k, dims = dims, 
                                  graph.name = graph_name, verbose = FALSE)
    
    graph <- Cell_cluster@graphs[[graph_name]]
    
    if (is.null(graph)) {
      stop("未生成图对象：", graph_name)
    }
    
    graph <- as(graph, "dgCMatrix")
    mean_connectivity <- mean(rowSums(graph > 0))
    
    results <- rbind(results, data.frame(k.param = k, mean_connectivity = mean_connectivity))
  }
  
  best_k <- results[which.max(results$mean_connectivity), "k.param"]
  cat("✅ best k.param:", best_k, "\n")
  
  return(list(best_k = best_k, results = results))
}


optimal_k <- find_optimal_k(Cell_cluster)
optimal_kparm <- optimal_k$best_k
Cell_cluster <- FindNeighbors(Cell_cluster, k.param = optimal_kparm, dims = 1:optimal_dim)

#三阶段分辨率选择
#TSNE选择
Cell_cluster <- RunTSNE(Cell_cluster, dims = 1:optimal_dim)
Cell_cluster <- RunUMAP(Cell_cluster, dims = 1:optimal_dim)
# 分辨率范围
resolutions <- seq(0.1, 1.6, by = 0.1)

# 创建数据框
metrics <- data.frame(
  resolution = resolutions,
  cluster_count = numeric(length(resolutions)),
  silhouette_score = numeric(length(resolutions)),
  marker_count = numeric(length(resolutions)),
  average_myAUC_all_clusters = numeric(length(resolutions))
)

# 计算函数，应用于每个分辨率
calculate_metrics <- function(res) {
  Cell_cluster <- FindClusters(Cell_cluster, resolution = res, verbose = FALSE)
  
  # 聚类数量
  cluster_count <- length(unique(Cell_cluster$seurat_clusters))
  
  # 计算轮廓系数
  dist_matrix <- dist(Cell_cluster@reductions$umap@cell.embeddings)
  sil <- silhouette(as.numeric(Cell_cluster$seurat_clusters), dist_matrix)
  silhouette_score <- mean(sil[, 3])
  
  # 标记基因数量
  markers <- FindAllMarkers(Cell_cluster, only.pos = TRUE, min.pct = 0.1, 
                            logfc.threshold = 0.5,
                            test.use = "roc")
  marker_count <- nrow(markers)
  
  # 计算前50个标记基因的平均 myAUC
  top_50_markers <- markers %>%
    group_by(cluster) %>%
    arrange(desc(myAUC)) %>%
    slice_head(n = 100) %>%
    ungroup()
  average_myAUC_all_clusters <- mean(top_50_markers$myAUC)
  
  # 返回结果
  return(c(cluster_count, silhouette_score, marker_count, average_myAUC_all_clusters))
}

# 使用 mclapply 并行计算每个分辨率的指标
results <- mclapply(resolutions, calculate_metrics, mc.cores = 2)


# 将结果存储到 metrics 数据框中
for (i in 1:length(resolutions)) {
  metrics[i, 2:5] <- results[[i]]
}

metrics$silhouette_score <- as.numeric(metrics$silhouette_score)
metrics$average_myAUC_all_clusters <- as.numeric(metrics$average_myAUC_all_clusters)
metrics$combined_score <- scale(metrics$silhouette_score) + scale(metrics$average_myAUC_all_clusters)

silhouette_scores <- metrics$silhouette_score
average_myAUC_all_clusters <- metrics$average_myAUC_all_clusters
combined_score <- metrics$combined_score
# 自定义梯形法函数
trapz <- function(x, y) {
  n <- length(x)
  area <- sum((x[2:n] - x[1:(n-1)]) * (y[2:n] + y[1:(n-1)]) / 2)
  return(area)
}

# 计算面积
area <- trapz(resolutions, silhouette_scores)
area_AUC <- trapz(resolutions, average_myAUC_all_clusters)
area_combined <- trapz(resolutions, combined_score)

# 梯度变化计算
gradients <- diff(silhouette_scores) / diff(resolutions)
gradients_AUC <- diff(average_myAUC_all_clusters) / diff(resolutions)
gradients_combined <- diff(combined_score) / diff(resolutions)
smooth_gradients <- stats::filter(gradients, rep(1/3, 3), sides = 2)
smooth_gradients_AUC <- stats::filter(gradients_AUC, rep(1/3, 3), sides = 2)
smooth_gradients_combined <- stats::filter(gradients_combined, rep(1/3, 3), sides = 2)
second_derivative <- diff(smooth_gradients, differences = 2)
second_derivative_AUC <- diff(smooth_gradients_AUC, differences = 2)
second_derivative_combined <- diff(smooth_gradients_combined, differences = 2)
elbow_point <- which.min(second_derivative) + 2
elbow_point_AUC <- which.min(second_derivative_AUC) + 2
elbow_point_combined <- which.min(second_derivative_combined) + 2
best_resolution <- resolutions[elbow_point]
best_resolution_AUC <- resolutions[elbow_point_AUC]
best_resolution_combined <- resolutions[elbow_point_combined]

pdf("silhouette_scores_plot.pdf", width = 8, height = 6)
plot(resolutions, silhouette_scores, type = "b", pch = 16, col = "blue",
     main = "Silhouette Scores with Resolutions",
     xlab = "Resolution", ylab = "Silhouette Score")
abline(v = best_resolution, col = "red", lty = 2)
text(best_resolution, silhouette_scores[elbow_point], labels = round(best_resolution, 2), pos = 4, col = "red")
dev.off()

pdf("myAUC_plot.pdf", width = 8, height = 6)
plot(resolutions, average_myAUC_all_clusters, type = "b", pch = 16, col = "blue",
     main = "Cluster myAUC with Resolutions",
     xlab = "Resolution", ylab = "Cluster myAUC")
abline(v = best_resolution_AUC, col = "red", lty = 2)
text(best_resolution_AUC, average_myAUC_all_clusters[elbow_point_AUC], labels = round(best_resolution_AUC, 2), pos = 4, col = "red")
dev.off()


pdf("combined_score_plot.pdf", width = 8, height = 6)
plot(resolutions, combined_score, type = "b", pch = 16, col = "blue",
     main = "combined_score with Resolutions",
     xlab = "Resolution", ylab = "Cluster combined_score")
abline(v = best_resolution_combined, col = "red", lty = 2)
text(best_resolution_combined, combined_score[elbow_point_combined], labels = round(best_resolution_combined, 2), pos = 4, col = "red")
dev.off()


clusters_A <- metrics$cluster_count[metrics$resolution == best_resolution]
clusters_B <- metrics$cluster_count[metrics$resolution == best_resolution_AUC]


if (best_resolution == best_resolution_AUC) {
  final_resolution <- best_resolution
} else {
  if (clusters_A == clusters_B) {
    final_resolution <- max(best_resolution, best_resolution_AUC)
  } else {
    final_resolution <- best_resolution_combined
  }
}

cat("final_resolution:", final_resolution, "\n")
Cell_Clusters <- FindClusters(Cell_cluster, resolution = final_resolution, verbose = FALSE)
rm(Cell_cluster)
save.image(file = "Resolution.RData")



