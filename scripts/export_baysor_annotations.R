# The purpose of this script is to export baysor segmentatation annotations to
# CSVs from the Seurat objects such that the annotations can be read in python
# Additionally, compute markers for each cell type and save as a CSV for access
# in python.

library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

obj <- readRDS(args[1])

md <- slot(object = obj, name = "meta.data")
md$cell <- rownames(md)
write.csv(md, file = args[2], row.names = FALSE)

Idents(obj) <- "predicted.class"
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.0, logfc.threshold = 0.0)
write.csv(markers, file = args[3], row.names = FALSE)
