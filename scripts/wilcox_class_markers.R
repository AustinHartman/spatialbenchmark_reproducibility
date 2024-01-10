# This script identifies markers using the linnarsson reference and save them to csv

library(Seurat)
library(ggplot2)
library(SeuratWrappers)

args <- commandArgs(trailingOnly = TRUE)

# Load reference and identify markers
ref <- readRDS(args[1])
ref <- NormalizeData(ref)
Idents(ref) <- "Class"
ref.markers <- RunPrestoAll(ref)
write.csv(ref.markers, args[2])

Idents(ref) <- "Taxonomy_group"
ref.markers <- RunPrestoAll(ref)
write.csv(ref.markers, args[3])
