# Save a heatmap

library(Seurat)
library(ggplot2)
library(SeuratWrappers)
source("./R/visualization.R")

args <- commandArgs(trailingOnly = TRUE)
message(paste("Arguments: ", args))

if (file.exists(args[1])) {
    obj <- readRDS(args[1])
} else {
    obj <- readRDS(file.path(getwd(), args[1]))
}
Idents(obj) <- "predicted.subclass"

markers <- RunPrestoAll(obj, slot = "data", only.pos = TRUE)
plot <- marker.heatmap(obj, markers, "predicted.subclass")

save.path <- file.path("./data/", paste0("heatmap_", args[2], ".png"))
message(paste("Saving to ", save.path))
ggsave(filename = save.path, plot = plot, units = "px", width = 3400, height = 3000)
