setwd("/DATA/")
load("STEP3.RData")
source("scoption.cfg")
library(Seurat)
library(clusterProfiler)
library(ggplot2)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(DESeq2)
library(KEGG.db)

dir.create("GOKEGG")
allcell <- levels(Cell_Clusters)
Cell_Clusters$celltype.Group <- paste(Idents(Cell_Clusters), Cell_Clusters$Group, sep = "_")
Idents(Cell_Clusters) <- "celltype.Group"
cell_list <- list()
for(i in 1:length(allcell)){
  current_celltype <- allcell[i]
  cell_list[[current_celltype]] <- Cell_Clusters[,Cell_Clusters$celltype == current_celltype, ]
  rm(current_celltype)
}

Comparison_list <- list()
DEG_list <- list()
for(i in 1:length(allcell)){
  current_celltype <- allcell[i]
  cell_group <- unique(cell_list[[current_celltype]]$celltype.Group)
  Comparison_list[[current_celltype]] <- paste0(cell_group[2],"_vs_",cell_group[1])
  DEG <- FindMarkers(cell_list[[current_celltype]], ident.1 = cell_group[1], ident.2 = cell_group[2], verbose = FALSE)
  DEG_up <- DEG[DEG$avg_log2FC > 0, ]
  DEG_down <- DEG[DEG$avg_log2FC < 0, ]
  DEG_list[[paste0(current_celltype, "_up")]] <- DEG_up
  write.csv(DEG_up,file.path("GOKEGG/" ,paste0(current_celltype, "_DEG_up.csv")))
  DEG_list[[paste0(current_celltype, "_down")]] <- DEG_down
  write.csv(DEG_down,file.path("GOKEGG/" ,paste0(current_celltype, "_DEG_down.csv")))
  rm(current_celltype, cell_group, DEG, DEG_up, DEG_down)
}
 


enrich_list <- list()
for(i in 1:length(allcell)){
  current_celltype <- allcell[i]
  enrich_list[[paste0(current_celltype, "_up")]] <- rownames(DEG_list[[paste0(current_celltype, "_up")]])
  enrich_list[[paste0(current_celltype, "_down")]] <- rownames(DEG_list[[paste0(current_celltype, "_down")]])
  rm(current_celltype)
}

enrich_entrezid_list <- list()
if (genome == "mouse"){
  for(i in 1:length(allcell)){
  current_celltype <- allcell[i]
  up=bitr(enrich_list[[paste0(current_celltype, "_up")]],'SYMBOL','ENTREZID','org.Mm.eg.db')
  down=bitr(enrich_list[[paste0(current_celltype, "_down")]],'SYMBOL','ENTREZID','org.Mm.eg.db')
  enrich_entrezid_list[[paste0(current_celltype, "_up")]] <- up$ENTREZID
  enrich_entrezid_list[[paste0(current_celltype, "_down")]] <- down$ENTREZID
  rm(current_celltype,up,down)
  }
}else if (genome == "human"){
  for(i in 1:length(allcell)){
    current_celltype <- allcell[i]
    up=bitr(enrich_list[[paste0(current_celltype, "_up")]],'SYMBOL','ENTREZID','org.Hs.eg.db')
    down=bitr(enrich_list[[paste0(current_celltype, "_down")]],'SYMBOL','ENTREZID','org.Hs.eg.db')
    enrich_entrezid_list[[paste0(current_celltype, "_up")]] <- up$ENTREZID
    enrich_entrezid_list[[paste0(current_celltype, "_down")]] <- down$ENTREZID
    rm(current_celltype,up,down)
  }
}else {
  system("echo 'Genome information not set'", intern = TRUE)
  stop("Please verify and set the genome. Exiting R.")
}

###GO富集
if (genome == "mouse"){
  for(i in 1:length(allcell)){
    current_celltype <- allcell[i]
    tryCatch({
    go_enrich_BP_up <- enrichGO(gene = enrich_entrezid_list[[paste0(current_celltype, "_up")]], OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP")
    go_enrich_BP_up_barplot <- barplot(go_enrich_BP_up , showCategory = 20, title = "GO Enrichment Analysis")
    go_enrich_BP_up_dotplot <- dotplot(go_enrich_BP_up , showCategory = 20, title = "GO Enrichment Analysis")
    ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_BP_up_barplot.pdf")),plot = go_enrich_BP_up_barplot,width = 7,height = 9)
    ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_BP_up_dotplot.pdf")),plot = go_enrich_BP_up_dotplot,width = 7,height = 9)
    write.table(go_enrich_BP_up, file.path("GOKEGG/" ,paste0(current_celltype, "_BP_up.txt")), sep = '\t', quote = FALSE, row.names = FALSE)
    go_enrich_BP_down <- enrichGO(gene = enrich_entrezid_list[[paste0(current_celltype, "_down")]], OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP")
    go_enrich_BP_down_barplot <- barplot(go_enrich_BP, showCategory = 20, title = "GO Enrichment Analysis")
    go_enrich_BP_down_dotplot <- dotplot(go_enrich_BP, showCategory = 20, title = "GO Enrichment Analysis")
    ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_BP_down_barplot.pdf")),plot = go_enrich_BP_down_barplot,width = 7,height = 9)
    ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_BP_down_dotplot.pdf")),plot = go_enrich_BP_down_dotplot,width = 7,height = 9)
    write.table(go_enrich_BP_down, file.path("GOKEGG/" ,paste0(current_celltype, "_BP_down.txt")), sep = '\t', quote = FALSE, row.names = FALSE)
    rm(go_enrich_BP_up,go_enrich_BP_up_barplot,go_enrich_BP_up_dotplot,go_enrich_BP_down,go_enrich_BP_down_barplot,go_enrich_BP_down_dotplot)
    go_enrich_CC_up <- enrichGO(gene = enrich_entrezid_list[[paste0(current_celltype, "_up")]], OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "CC")
    go_enrich_CC_up_barplot <- barplot(go_enrich_CC_up , showCategory = 20, title = "GO Enrichment Analysis")
    go_enrich_CC_up_dotplot <- dotplot(go_enrich_CC_up , showCategory = 20, title = "GO Enrichment Analysis")
    ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_CC_up_barplot.pdf")),plot = go_enrich_CC_up_barplot,width = 7,height = 9)
    ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_CC_up_dotplot.pdf")),plot = go_enrich_CC_up_dotplot,width = 7,height = 9)
    write.table(go_enrich_CC_up, file.path("GOKEGG/" ,paste0(current_celltype, "_CC_up.txt")), sep = '\t', quote = FALSE, row.names = FALSE)
    go_enrich_CC_down <- enrichGO(gene = enrich_entrezid_list[[paste0(current_celltype, "_down")]], OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "CC")
    go_enrich_CC_down_barplot <- barplot(go_enrich_CC, showCategory = 20, title = "GO Enrichment Analysis")
    go_enrich_CC_down_dotplot <- dotplot(go_enrich_CC, showCategory = 20, title = "GO Enrichment Analysis")
    ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_CC_down_barplot.pdf")),plot = go_enrich_CC_down_barplot,width = 7,height = 9)
    ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_CC_down_dotplot.pdf")),plot = go_enrich_CC_down_dotplot,width = 7,height = 9)
    write.table(go_enrich_CC_down, file.path("GOKEGG/" ,paste0(current_celltype, "_CC_down.txt")), sep = '\t', quote = FALSE, row.names = FALSE)
    rm(go_enrich_CC_up,go_enrich_CC_up_barplot,go_enrich_CC_up_dotplot,go_enrich_CC_down,go_enrich_CC_down_barplot,go_enrich_CC_down_dotplot)
    go_enrich_MF_up <- enrichGO(gene = enrich_entrezid_list[[paste0(current_celltype, "_up")]], OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "MF")
    go_enrich_MF_up_barplot <- barplot(go_enrich_MF_up , showCategory = 20, title = "GO Enrichment Analysis")
    go_enrich_MF_up_dotplot <- dotplot(go_enrich_MF_up , showCategory = 20, title = "GO Enrichment Analysis")
    ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_MF_up_barplot.pdf")),plot = go_enrich_MF_up_barplot,width = 7,height = 9)
    ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_MF_up_dotplot.pdf")),plot = go_enrich_MF_up_dotplot,width = 7,height = 9)
    write.table(go_enrich_MF_up, file.path("GOKEGG/" ,paste0(current_celltype, "_MF_up.txt")), sep = '\t', quote = FALSE, row.names = FALSE)
    go_enrich_MF_down <- enrichGO(gene = enrich_entrezid_list[[paste0(current_celltype, "_down")]], OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "MF")
    go_enrich_MF_down_barplot <- barplot(go_enrich_MF, showCategory = 20, title = "GO Enrichment Analysis")
    go_enrich_MF_down_dotplot <- dotplot(go_enrich_MF, showCategory = 20, title = "GO Enrichment Analysis")
    ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_MF_down_barplot.pdf")),plot = go_enrich_MF_down_barplot,width = 7,height = 9)
    ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_MF_down_dotplot.pdf")),plot = go_enrich_MF_down_dotplot,width = 7,height = 9)
    write.table(go_enrich_MF_down, file.path("GOKEGG/" ,paste0(current_celltype, "_MF_down.txt")), sep = '\t', quote = FALSE, row.names = FALSE)
    rm(go_enrich_MF_up,go_enrich_MF_up_barplot,go_enrich_MF_up_dotplot,go_enrich_MF_down,go_enrich_MF_down_barplot,go_enrich_MF_down_dotplot)
    }, error = function(e) {
      cat("Error occurred. Skipping this iteration.\n")
    })
     }
}else if (genome == "human"){
  for(i in 1:length(allcell)){
    current_celltype <- allcell[i]
    tryCatch({
      go_enrich_BP_up <- enrichGO(gene = enrich_entrezid_list[[paste0(current_celltype, "_up")]], OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = "BP")
      go_enrich_BP_up_barplot <- barplot(go_enrich_BP_up , showCategory = 20, title = "GO Enrichment Analysis")
      go_enrich_BP_up_dotplot <- dotplot(go_enrich_BP_up , showCategory = 20, title = "GO Enrichment Analysis")
      ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_BP_up_barplot.pdf")),plot = go_enrich_BP_up_barplot,width = 7,height = 9)
      ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_BP_up_dotplot.pdf")),plot = go_enrich_BP_up_dotplot,width = 7,height = 9)
      write.table(go_enrich_BP_up, file.path("GOKEGG/" ,paste0(current_celltype, "_BP_up.txt")), sep = '\t', quote = FALSE, row.names = FALSE)
      go_enrich_BP_down <- enrichGO(gene = enrich_entrezid_list[[paste0(current_celltype, "_down")]], OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = "BP")
      go_enrich_BP_down_barplot <- barplot(go_enrich_BP, showCategory = 20, title = "GO Enrichment Analysis")
      go_enrich_BP_down_dotplot <- dotplot(go_enrich_BP, showCategory = 20, title = "GO Enrichment Analysis")
      ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_BP_down_barplot.pdf")),plot = go_enrich_BP_down_barplot,width = 7,height = 9)
      ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_BP_down_dotplot.pdf")),plot = go_enrich_BP_down_dotplot,width = 7,height = 9)
      write.table(go_enrich_BP_down, file.path("GOKEGG/" ,paste0(current_celltype, "_BP_down.txt")), sep = '\t', quote = FALSE, row.names = FALSE)
      rm(go_enrich_BP_up,go_enrich_BP_up_barplot,go_enrich_BP_up_dotplot,go_enrich_BP_down,go_enrich_BP_down_barplot,go_enrich_BP_down_dotplot)
      go_enrich_CC_up <- enrichGO(gene = enrich_entrezid_list[[paste0(current_celltype, "_up")]], OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = "CC")
      go_enrich_CC_up_barplot <- barplot(go_enrich_CC_up , showCategory = 20, title = "GO Enrichment Analysis")
      go_enrich_CC_up_dotplot <- dotplot(go_enrich_CC_up , showCategory = 20, title = "GO Enrichment Analysis")
      ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_CC_up_barplot.pdf")),plot = go_enrich_CC_up_barplot,width = 7,height = 9)
      ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_CC_up_dotplot.pdf")),plot = go_enrich_CC_up_dotplot,width = 7,height = 9)
      write.table(go_enrich_CC_up, file.path("GOKEGG/" ,paste0(current_celltype, "_CC_up.txt")), sep = '\t', quote = FALSE, row.names = FALSE)
      go_enrich_CC_down <- enrichGO(gene = enrich_entrezid_list[[paste0(current_celltype, "_down")]], OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = "CC")
      go_enrich_CC_down_barplot <- barplot(go_enrich_CC, showCategory = 20, title = "GO Enrichment Analysis")
      go_enrich_CC_down_dotplot <- dotplot(go_enrich_CC, showCategory = 20, title = "GO Enrichment Analysis")
      ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_CC_down_barplot.pdf")),plot = go_enrich_CC_down_barplot,width = 7,height = 9)
      ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_CC_down_dotplot.pdf")),plot = go_enrich_CC_down_dotplot,width = 7,height = 9)
      write.table(go_enrich_CC_down, file.path("GOKEGG/" ,paste0(current_celltype, "_CC_down.txt")), sep = '\t', quote = FALSE, row.names = FALSE)
      rm(go_enrich_CC_up,go_enrich_CC_up_barplot,go_enrich_CC_up_dotplot,go_enrich_CC_down,go_enrich_CC_down_barplot,go_enrich_CC_down_dotplot)
      go_enrich_MF_up <- enrichGO(gene = enrich_entrezid_list[[paste0(current_celltype, "_up")]], OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = "MF")
      go_enrich_MF_up_barplot <- barplot(go_enrich_MF_up , showCategory = 20, title = "GO Enrichment Analysis")
      go_enrich_MF_up_dotplot <- dotplot(go_enrich_MF_up , showCategory = 20, title = "GO Enrichment Analysis")
      ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_MF_up_barplot.pdf")),plot = go_enrich_MF_up_barplot,width = 7,height = 9)
      ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_MF_up_dotplot.pdf")),plot = go_enrich_MF_up_dotplot,width = 7,height = 9)
      write.table(go_enrich_MF_up, file.path("GOKEGG/" ,paste0(current_celltype, "_MF_up.txt")), sep = '\t', quote = FALSE, row.names = FALSE)
      go_enrich_MF_down <- enrichGO(gene = enrich_entrezid_list[[paste0(current_celltype, "_down")]], OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = "MF")
      go_enrich_MF_down_barplot <- barplot(go_enrich_MF, showCategory = 20, title = "GO Enrichment Analysis")
      go_enrich_MF_down_dotplot <- dotplot(go_enrich_MF, showCategory = 20, title = "GO Enrichment Analysis")
      ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_MF_down_barplot.pdf")),plot = go_enrich_MF_down_barplot,width = 7,height = 9)
      ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_MF_down_dotplot.pdf")),plot = go_enrich_MF_down_dotplot,width = 7,height = 9)
      write.table(go_enrich_MF_down, file.path("GOKEGG/" ,paste0(current_celltype, "_MF_down.txt")), sep = '\t', quote = FALSE, row.names = FALSE)
      rm(go_enrich_MF_up,go_enrich_MF_up_barplot,go_enrich_MF_up_dotplot,go_enrich_MF_down,go_enrich_MF_down_barplot,go_enrich_MF_down_dotplot)
    }, error = function(e) {
      cat("Error occurred. Skipping this iteration.\n")
    })
    }
}else {
  system("echo 'Genome information not set'", intern = TRUE)
  stop("Please verify and set the genome. Exiting R.")
}



if (genome == "mouse"){
  for(i in 1:length(allcell)){
    current_celltype <- allcell[i]
    tryCatch({
    kegg_enrich_up <- enrichKEGG(gene = enrich_entrezid_list[[paste0(current_celltype, "_up")]], organism = 'mmu', keyType = 'kegg', pvalueCutoff = 1,use_internal_data = TRUE)
    enrich_KEGG_up_barplot <- barplot(kegg_enrich_up , showCategory = 20, title = "KEGG Pathway")
    enrich_KEGG_up_dotplot <- dotplot(kegg_enrich_up , showCategory = 20, title = "KEGG Pathway")
    ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_KEGG_up_barplot.pdf")),plot = enrich_KEGG_up_barplot,width = 7,height = 9)
    ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_KEGG_up_dotplot.pdf")),plot = enrich_KEGG_up_barplot,width = 7,height = 9)
    write.table(kegg_enrich_up, file.path("GOKEGG/" ,paste0(current_celltype, "_KEGG_up.txt")), sep = '\t', quote = FALSE, row.names = FALSE)
    kegg_enrich_down <- enrichKEGG(gene = enrich_entrezid_list[[paste0(current_celltype, "_down")]], organism = 'mmu', keyType = 'kegg', pvalueCutoff = 1,use_internal_data = TRUE)
    enrich_KEGG_down_barplot <- barplot(kegg_enrich_down , showCategory = 20, title = "KEGG Pathway")
    enrich_KEGG_down_dotplot <- dotplot(kegg_enrich_down , showCategory = 20, title = "KEGG Pathway")
    ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_KEGG_down_barplot.pdf")),plot = enrich_KEGG_down_barplot,width = 7,height = 9)
    ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_KEGG_down_dotplot.pdf")),plot = enrich_KEGG_down_barplot,width = 7,height = 9)
    write.table(kegg_enrich_down, file.path("GOKEGG/" ,paste0(current_celltype, "_KEGG_down.txt")), sep = '\t', quote = FALSE, row.names = FALSE)
    rm(kegg_enrich_up,enrich_KEGG_up_barplot,enrich_KEGG_up_dotplot,kegg_enrich_down,enrich_KEGG_down_barplot,enrich_KEGG_down_barplot)
    }, error = function(e) {
      cat("Error occurred. Skipping this iteration.\n")
    })
  }
}else if (genome == "human"){
  for(i in 1:length(allcell)){
    current_celltype <- allcell[i]
    tryCatch({
      kegg_enrich_up <- enrichKEGG(gene = enrich_entrezid_list[[paste0(current_celltype, "_up")]], organism = 'hsa', keyType = 'kegg', pvalueCutoff = 1,use_internal_data = TRUE)
      enrich_KEGG_up_barplot <- barplot(kegg_enrich_up , showCategory = 20, title = "KEGG Pathway")
      enrich_KEGG_up_dotplot <- dotplot(kegg_enrich_up , showCategory = 20, title = "KEGG Pathway")
      ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_KEGG_up_barplot.pdf")),plot = enrich_KEGG_up_barplot,width = 7,height = 9)
      ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_KEGG_up_dotplot.pdf")),plot = enrich_KEGG_up_barplot,width = 7,height = 9)
      write.table(kegg_enrich_up, file.path("GOKEGG/" ,paste0(current_celltype, "_KEGG_up.txt")), sep = '\t', quote = FALSE, row.names = FALSE)
      kegg_enrich_down <- enrichKEGG(gene = enrich_entrezid_list[[paste0(current_celltype, "_down")]], organism = 'hsa', keyType = 'kegg', pvalueCutoff = 1,use_internal_data = TRUE)
      enrich_KEGG_down_barplot <- barplot(kegg_enrich_down , showCategory = 20, title = "KEGG Pathway")
      enrich_KEGG_down_dotplot <- dotplot(kegg_enrich_down , showCategory = 20, title = "KEGG Pathway")
      ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_KEGG_down_barplot.pdf")),plot = enrich_KEGG_down_barplot,width = 7,height = 9)
      ggsave(file.path("GOKEGG/" ,paste0(current_celltype, "_KEGG_down_dotplot.pdf")),plot = enrich_KEGG_down_barplot,width = 7,height = 9)
      write.table(kegg_enrich_down, file.path("GOKEGG/" ,paste0(current_celltype, "_KEGG_down.txt")), sep = '\t', quote = FALSE, row.names = FALSE)
      rm(kegg_enrich_up,enrich_KEGG_up_barplot,enrich_KEGG_up_dotplot,kegg_enrich_down,enrich_KEGG_down_barplot,enrich_KEGG_down_barplot)
    }, error = function(e) {
      cat("Error occurred. Skipping this iteration.\n")
    })
  }
}else {
  system("echo 'Genome information not set'", intern = TRUE)
  stop("Please verify and set the genome. Exiting R.")
}

save.image(file = "STEP4.RData")
