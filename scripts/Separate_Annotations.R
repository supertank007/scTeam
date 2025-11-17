library(Seurat)
library(R.utils)
library(data.table)

setwd(Sys.getenv("PWD"))
geneinfo<- readRDS('/home/code/geneinfo.rds')
load("Resolution.RData")
cluster_list <- SplitObject(Cell_cluster, split.by = "seurat_clusters")
cluster_rna_data <- list()

for (i in 1:length(cluster_list)) {
  rna_data <- cluster_list[[i]][['RNA']]@data
  cluster_rna_data[[paste0("cluster_", i)]] <- rna_data
  cat("Stored RNA data for cluster", i, "\n")
}



for (m in 1:length(cluster_list)) {
genename<- rownames(cluster_rna_data[[m]])
genename1<- genename[genename %in% geneinfo$Symbol]
genename2<- genename[!genename %in% geneinfo$Symbol]
genename3<- genename2[genename2 %in% geneinfo$Synonyms]
genename4<- rep('NA',length(genename3))
for (i in 1:length(genename3)) {
  d1<- geneinfo[geneinfo$Synonyms == genename3[i],]$Symbol
  if(length(d1) == 1){
    genename4[i]<- d1
  }
}
genename3<- c(genename1,genename3)
genename4<- c(genename1,genename4)
genedata<- data.frame(raw_name = genename3,new_name = genename4,stringsAsFactors = F)
genedata<- genedata[!genedata$new_name == 'NA',]
genedata1<- as.data.frame(table(genedata$new_name),stringsAsFactors = F)
genedata1<- genedata1[genedata1$Freq == 1,]
genedata<- genedata[genedata$new_name %in% genedata1$Var1,]

cluster_rna_data[[m]] <- cluster_rna_data[[m]][genedata$raw_name,]
all(rownames(cluster_rna_data[[m]]) == genedata$raw_name)
rownames(cluster_rna_data[[m]])<- genedata$new_name
all(rownames(cluster_rna_data[[m]]) == genedata$new_name)
all(rownames(cluster_rna_data[[m]]) %in% geneinfo$Symbol)

cluster_rna_data[[m]]<- as.matrix(cluster_rna_data[[m]])
file_name <- paste0("cluster_", m, ".csv")
write.csv(cluster_rna_data[[m]], file = file_name)
cat("Saved RNA data for cluster", m, "as", file_name, "\n")

}

