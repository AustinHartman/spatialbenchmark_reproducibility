# Purpose of this script is to convert high confidence matrices to Seurat objects

library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3) {
    cutoff <- as.numeric(args[4])
} else {
   cutoff <- 0.0
}

df <- read.csv(args[1], header = TRUE, stringsAsFactors = FALSE)
df <- df[nchar(df$cell) > 0, ]

# lowering the parameter increases sensitivity with a decrease in specificity
df <- df[df$assignment_confidence >= cutoff, ]
mtx <- table(df$cell, df$gene)
mtx <- as(mtx, Class = "Matrix")

library(Seurat)
obj <- CreateSeuratObject(t(mtx))

save.path <- file.path("data", paste0(args[2], "_baysor_", args[3], ".rds"))
saveRDS(obj, file = save.path)
