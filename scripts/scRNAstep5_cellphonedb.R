setwd("/DATA/")
load("STEP3.RData")
library(Seurat)
library(data.table)
library(AnnoProbe) 
source("scoption.cfg")

if (genome == "mouse") {
  exp  <- as.data.frame(Cell_Clusters@assays$RNA@data)
  rownames(exp) <- toupper(rownames(exp))
  ids=annoGene( rownames(exp),'SYMBOL','human')
  length(unique(ids$ENSEMBL)) 
  ids=ids[!duplicated(ids$ENSEMBL),]
  length(unique(ids$SYMBOL)) 
  ids=ids[!duplicated(ids$SYMBOL),]
  cpdb_counts=exp[ids$SYMBOL,]
  cpdb_counts[1:4,1:4]
  table(Cell_Clusters$celltype ) 
  cpdb_meta <-  data.frame(Cell = rownames(Cell_Clusters@meta.data), 
                           cell_type = Cell_Clusters$celltype )  
  head(cpdb_meta)
  cpdb_meta$cell_type=gsub(' ','_',cpdb_meta$cell_type)
  cpdb_meta$cell_type=gsub('\\+','',cpdb_meta$cell_type) 
  table(cpdb_meta$cell_type)
  length(unique(cpdb_meta$Cell))
  identical(colnames(cpdb_counts),cpdb_meta$Cell)
  rownames(cpdb_counts)=ids$ENSEMBL
  cpdb_counts=cbind(rownames(cpdb_counts),cpdb_counts)
  colnames(cpdb_counts)[1]='Gene'
  cpdb_counts[1:4,1:4]
  
  fwrite(cpdb_counts, file = "cpdb_counts.txt", quote = FALSE, sep = "\t", row.names = FALSE)
  fwrite(cpdb_meta, file = "cpdb_meta.txt", row.names=FALSE, sep='\t',quote = FALSE) 
} else if (genome == "human") {
  exp  <- as.data.frame(Cell_Clusters@assays$RNA@data)
  ids=annoGene( rownames(exp),'SYMBOL','human')
  length(unique(ids$ENSEMBL)) 
  ids=ids[!duplicated(ids$ENSEMBL),]
  length(unique(ids$SYMBOL)) 
  ids=ids[!duplicated(ids$SYMBOL),]
  cpdb_counts=exp[ids$SYMBOL,]
  cpdb_counts[1:4,1:4]
  table(Cell_Clusters$celltype ) 
  cpdb_meta <-  data.frame(Cell = rownames(Cell_Clusters@meta.data), 
                           cell_type = Cell_Clusters$celltype )  
  head(cpdb_meta)
  cpdb_meta$cell_type=gsub(' ','_',cpdb_meta$cell_type)
  cpdb_meta$cell_type=gsub('\\+','',cpdb_meta$cell_type) 
  table(cpdb_meta$cell_type)
  length(unique(cpdb_meta$Cell))
  identical(colnames(cpdb_counts),cpdb_meta$Cell)
  rownames(cpdb_counts)=ids$ENSEMBL
  cpdb_counts=cbind(rownames(cpdb_counts),cpdb_counts)
  colnames(cpdb_counts)[1]='Gene'
  cpdb_counts[1:4,1:4]
  
  fwrite(cpdb_counts, file = "cpdb_counts.txt", quote = FALSE, sep = "\t", row.names = FALSE)
  fwrite(cpdb_meta, file = "cpdb_meta.txt", row.names=FALSE, sep='\t',quote = FALSE)
} else {
  system("echo 'Genome information not set'", intern = TRUE)
  stop("Please verify and set the genome. Exiting R.")
}
