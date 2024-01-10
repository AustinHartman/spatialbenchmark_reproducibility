# Purpose of this script is to convert 10x Xenium molecule files to Seurat objects

library(Matrix)
args <- commandArgs(trailingOnly = TRUE)

df <- read.csv(args[1], header = TRUE, stringsAsFactors = FALSE)
df <- df[df$cell > 0, ]
mtx <- table(df$cell, df$gene)
mtx <- as(mtx, Class = "Matrix")

library(Seurat)
obj <- CreateSeuratObject(t(mtx))

save.path <- file.path("data", paste0(args[2], "_cortex.rds"))
saveRDS(obj, file = save.path)
