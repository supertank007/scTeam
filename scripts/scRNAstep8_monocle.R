setwd("/DATA/")
load("STEP3.RData")
library(Seurat)
library(monocle)
library(ggpubr)
library(patchwork)
library(dplyr)
library(ggplot2)
source("scoption.cfg")

dir.create("monocle")
allcell <- levels(Cell_Clusters)
cell_list <- list()
for(i in 1:length(allcell)){
  current_celltype <- allcell[i]
  cell_list[[current_celltype]] <- Cell_Clusters[,Cell_Clusters$celltype == current_celltype, ]
  rm(current_celltype)
}


for(i in 1:length(allcell)){
  current_celltype <- allcell[i]
  tryCatch({
  data <- as(as.matrix(cell_list[[current_celltype]]@assays$RNA@counts), 'sparseMatrix')
  pd <- new('AnnotatedDataFrame', data = cell_list[[current_celltype]]@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table,mean_expression >= 0.1 & dispersion_empirical >= 1 *dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
order_genes <- plot_ordering_genes(cds)
ggsave(file.path("monocle/", paste0(current_celltype,"_order_genes.pdf")),plot = order_genes)

cds <- reduceDimension(
  cds,
  max_components = 2,
  method = 'DDRTree')

cds <- orderCells(cds)

colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")


a1 <- plot_cell_trajectory(cds, color_by = "celltype") + scale_color_manual(values = colour)

a2 <- plot_cell_trajectory(cds, color_by = "State") + scale_color_manual(values = colour)

a3 <- plot_cell_trajectory(cds, color_by = "Pseudotime") 

ggsave(file.path("monocle/", paste0(current_celltype,"_celltype.pdf")),plot = a1)
ggsave(file.path("monocle/", paste0(current_celltype,"_State.pdf")),plot = a2)
ggsave(file.path("monocle/", paste0(current_celltype,"_Pseudotime.pdf")),plot = a3)

Pseudotime_tree <- plot_complex_cell_trajectory(cds, x = 1, y = 2,color_by = "Pseudotime")
ggsave(file.path("monocle/", paste0(current_celltype,"_Pseudotime_tree.pdf")),plot = Pseudotime_tree,width = 9,height = 6)


df <- pData(cds)

Pseudotime_cell <- ggplot(df, aes(Pseudotime, colour = State, fill=State)) +
  
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()

ggsave(file.path("monocle/", paste0(current_celltype,"_Pseudotime_cell.pdf")),plot = Pseudotime_cell,width = 8,height = 5)

keygenes <- head(disp.genes,5)

cds_subset <- cds[keygenes,]

p1 <-plot_genes_in_pseudotime(cds_subset, color_by = "State")

p2 <-plot_genes_in_pseudotime(cds_subset, color_by = "celltype")

p3 <-plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")

ggsave(file.path("monocle/", paste0(current_celltype,"_State_genetop5.pdf")),plot = p1,width = 5,height = 10)
ggsave(file.path("monocle/", paste0(current_celltype,"_celltype_genetop5.pdf")),plot = p2,width = 5,height = 10)
ggsave(file.path("monocle/", paste0(current_celltype,"_Pseudotime_genetop5.pdf")),plot = p3,width = 5,height = 10)

#############

if (current_celltype %in% names(Monocle_list)) {
  for(i in 1:length(Monocle_list[[current_celltype]])){
    current_gene <- Monocle_list[[current_celltype]][i]
    pData(cds)[, current_gene] = log2(exprs(cds)[current_gene, ] + 1)
    p <- plot_cell_trajectory(cds, color_by = current_gene)
    ggsave(file.path("monocle/", paste0(current_celltype,"_Pseudotime_", current_gene, ".pdf")),plot = p,width = 4,height = 4)
  cat("Found", current_gene,"in Monocle_list")}
} else {
  cat("Table", current_celltype, "not found in Monocle_list")
}

########################

Time_diff <-differentialGeneTest(cds[disp.genes,], cores = 1,
                                 
                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")

Time_diff_top50 <- Time_diff[order(Time_diff$pval), ]
Time_diff_top50 <- head(Time_diff_top50, 50)

write.csv(Time_diff,file.path("monocle/", paste0(current_celltype,"Time_diff_all.csv")), row.names = F)

Time_genes <- Time_diff %>%pull(gene_short_name) %>% as.character()
heat1 <- plot_pseudotime_heatmap(cds[Time_genes,], show_rownames=T, return_heatmap=T)
ggsave(file.path("monocle/", paste0(current_celltype,"_heatmap_all.pdf")),plot = heat1,width = 5,height = 8)
Time_genes <- Time_diff_top50 %>%pull(gene_short_name) %>% as.character()
heat2 <- plot_pseudotime_heatmap(cds[Time_genes,], show_rownames=T, return_heatmap=T)
ggsave(file.path("monocle/", paste0(current_celltype,"_heatmap_top50.pdf")),plot = heat2,width = 5,height = 8)

}, error = function(e) {
  cat("Error occurred. Skipping this iteration.\n")
})
}


###all in one
data <- as(as.matrix(Cell_Clusters@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = Cell_Clusters@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table,mean_expression >= 0.1 & dispersion_empirical >= 1 *dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
order_genes <- plot_ordering_genes(cds)
ggsave(file.path("monocle/", "all_order_genes.pdf"),plot = order_genes)
cds <- reduceDimension(
  cds,
  max_components = 2,
  method = 'DDRTree')
print("order all cell please wait.....This may take a lot of time")
cds <- orderCells(cds)
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
a1 <- plot_cell_trajectory(cds, color_by = "celltype") + scale_color_manual(values = colour)
a2 <- plot_cell_trajectory(cds, color_by = "State") + scale_color_manual(values = colour)
a3 <- plot_cell_trajectory(cds, color_by = "Pseudotime") 
ggsave(file.path("monocle/", "all_celltype.pdf"), plot = a1)
ggsave(file.path("monocle/", "all_State.pdf"), plot = a2)
ggsave(file.path("monocle/", "all_Pseudotime.pdf"), plot = a3)
Pseudotime_tree <- plot_complex_cell_trajectory(cds, x = 1, y = 2,color_by = "Pseudotime")
Pseudotime_tree_type <- plot_complex_cell_trajectory(cds, x = 1, y = 2,color_by = "celltype")
ggsave(file.path("monocle/", "all_Pseudotime_tree.pdf"), plot = Pseudotime_tree, width = 9,height = 6)
ggsave(file.path("monocle/", "all_Pseudotime_tree_type.pdf"), plot = Pseudotime_tree_type, width = 9,height = 6)
df <- pData(cds)
Pseudotime_cell <- ggplot(df, aes(Pseudotime, colour = celltype, fill=celltype)) + geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
ggsave(file.path("monocle/", "all_Pseudotime_cell.pdf"),plot = Pseudotime_cell,width = 8,height = 5)
keygenes <- head(disp.genes,5)
cds_subset <- cds[keygenes,]
p1 <-plot_genes_in_pseudotime(cds_subset, color_by = "State")
p2 <-plot_genes_in_pseudotime(cds_subset, color_by = "celltype")
p3 <-plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
ggsave(file.path("monocle/", "all_State_genetop5.pdf"),plot = p1,width = 5,height = 10)
ggsave(file.path("monocle/", "all_celltype_genetop5.pdf"),plot = p2,width = 5,height = 10)
ggsave(file.path("monocle/", "all_Pseudotime_genetop5.pdf"),plot = p3,width = 5,height = 10)
Time_diff <-differentialGeneTest(cds[disp.genes,], cores = 1,
                                 
                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff_top50 <- Time_diff[order(Time_diff$pval), ]
Time_diff_top50 <- head(Time_diff_top50, 50)
write.csv(Time_diff,file.path("monocle/", "_all_Time_diff_all.csv"), row.names = F)
Time_genes <- Time_diff %>%pull(gene_short_name) %>% as.character()
heat1 <- plot_pseudotime_heatmap(cds[Time_genes,], show_rownames=T, return_heatmap=T)
ggsave(file.path("monocle/", "all_heatmap_all.pdf"),plot = heat1,width = 5,height = 8)
Time_genes <- Time_diff_top50 %>%pull(gene_short_name) %>% as.character()
heat2 <- plot_pseudotime_heatmap(cds[Time_genes,], show_rownames=T, return_heatmap=T)
ggsave(file.path("monocle/", "all_heatmap_top50.pdf"),plot = heat2,width = 5,height = 8)
