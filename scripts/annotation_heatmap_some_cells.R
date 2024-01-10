# Save a heatmap

library(Seurat)
library(ggplot2)
library(SeuratWrappers)
source("./R/visualization.R")

args <- commandArgs(trailingOnly = TRUE)
message(paste("Arguments: ", args))
group.by <- "predicted.subclass"
if (args[2] == "linnarsson") group.by <- "Subclass"

marker.heatmap <- function(
    object,
    markers,
    group.by
) {
  Idents(object) <- group.by
  heatmap.genes <- c()
  for (ct in c("Neurons", "Oligos", "Astrocyte", "Vascular", "Immune")) {
    if (ct %in% Idents(object)) {
      heatmap.genes <- c(
        heatmap.genes,
        head(markers[markers$cluster == ct, ]$gene, 5)
      )
    }
  }
  cells <- c()
  for (ct in c("Neurons", "Oligos", "Astrocyte", "Vascular", "Immune")) {
    if (ct %in% Idents(object)) {
      ct.cells <- WhichCells(object, idents = ct)
      if (length(ct.cells) > 50) {
        cells <- c(cells, sample(ct.cells, 50))
      } else {
        cells <- c(cells, ct.cells)
      }
    }
  }
  object <- subset(object, cells = cells)
  Idents(object) <- group.by
  object <- ScaleData(object, features = heatmap.genes)
  object@assays$RNA@counts[object@assays$RNA@counts > 5] <- 5
  object.heatmap <- DoHeatmap(
    object, features = heatmap.genes, cells = cells, disp.min=0, disp.max = 5, slot = "counts"
  ) + ggtitle("Heatmap")
  return(object.heatmap)
}

if (file.exists(args[1])) {
    obj <- readRDS(args[1])
} else {
    obj <- readRDS(file.path(getwd(), args[1]))
}
Idents(obj) <- group.by

markers <- RunPrestoAll(obj, slot = "data", only.pos = TRUE)
plot <- marker.heatmap(obj, markers, group.by) +
  theme(text = element_text(size = 26, colour = "black")) +
  ggtitle("") + NoLegend()

save.path <- file.path("./data/", paste0("heatmap_subset_", args[2], ".png"))
message(paste("Saving to ", save.path))
ggsave(filename = save.path, plot = plot, units = "px", width = 3000, height = 3000)
