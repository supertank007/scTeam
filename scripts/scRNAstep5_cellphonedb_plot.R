setwd("/DATA/")
library(CellChat)
library(tidyr)
df.net <- read.table("out/count_network.txt",
                     header = T,sep = "\t",stringsAsFactors = F)
meta.data <- read.table("cpdb_meta.txt",
                        header = T,sep = "\t",stringsAsFactors = F)

groupSize <- as.numeric(table(meta.data$cell_type))

df.net <- spread(df.net, TARGET, count)
rownames(df.net) <- df.net$SOURCE
df.net <- df.net[, -1]
df.net <- as.matrix(df.net)

pdf(file.path("out/", "netVisual_circle.pdf"))
netVisual_circle(df.net, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
dev.off()

# par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(df.net)) {
  mat2 <- matrix(0, nrow = nrow(df.net), ncol = ncol(df.net), dimnames = dimnames(df.net))
  mat2[i, ] <- df.net[i, ]
  file_name <- paste0("netVisual_circle_", i, ".pdf")
  pdf(file.path("out/", file_name))  
  netVisual_circle(mat2, 
                   vertex.weight = groupSize,
                   weight.scale = T, 
                   edge.weight.max = max(df.net),
                   title.name = rownames(df.net)[i],
                   arrow.size=0.5)
  dev.off() 
}
