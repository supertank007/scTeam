setwd(Sys.getenv("PWD"))
load("STEP3.RData")
library(Seurat)
meta <- Cell_Clusters@meta.data
meta$barcode <- gsub("_\\d+$", "", colnames(Cell_Clusters))
barcode_list <- split(meta$barcode, meta$id)
output_dir <- "barcodes_by_sample"
dir.create(output_dir, showWarnings = FALSE)

for (sample in names(barcode_list)) {
  barcodes <- barcode_list[[sample]]
  file_path <- file.path(output_dir, paste0(sample, ".barcodes.txt"))
  writeLines(barcodes, file_path)
}
