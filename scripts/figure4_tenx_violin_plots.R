library(Seurat)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

tenx <- readRDS(args[1])
tenx.annotated <- readRDS(args[2])
tenx$annotations <- tenx.annotated$predicted.class

cell.names <- tenx@images$fov@boundaries$centroids@cells
pos <- tenx@images$fov@boundaries$centroids@coords
df <- as.data.frame(pos)
df$cell <- cell.names

thalamus.cells <- df[
  (df$x>3100 & df$x<4400 & df$y>2800 & df$y<4100)
  , ]
cortex.cells <- df[
  (df$x>3300 & df$x<4300 & df$y>5600 & df$y<6600)
  , ]

tenx$annotations <- ifelse(tenx$annotations == "Astrocytes" & Cells(tenx) %in% thalamus.cells$cell, "Astrocytes_thalamus", tenx$annotations)
tenx$annotations <- ifelse(tenx$annotations == "Astrocytes" & Cells(tenx) %in% cortex.cells$cell, "Astrocytes_cortex", tenx$annotations)
tenx$annotations <- ifelse(tenx$annotations == "Neurons" & Cells(tenx) %in% thalamus.cells$cell, "Neurons_thalamus", tenx$annotations)
tenx$annotations <- ifelse(tenx$annotations == "Neurons" & Cells(tenx) %in% cortex.cells$cell, "Neurons_cortex", tenx$annotations)

Idents(tenx) <- "annotations"
Idents(tenx) <- as.factor(tenx$annotations)
tenx <- NormalizeData(tenx)
p1 <- VlnPlot(
  tenx, idents=c("Astrocytes_cortex", "Astrocytes_thalamus", "Neurons_cortex", "Neurons_thalamus"),
  features=c("Slc17a6"), pt.size=0, ncol=1, cols = c("#d55e00", "#d55e00", "#d55e00", "#d55e00")
) + NoLegend() + theme(text=element_text(size=30))
p2 <- VlnPlot(
  tenx, idents=c("Astrocytes_cortex", "Astrocytes_thalamus", "Neurons_cortex", "Neurons_thalamus"),
  features=c("Satb2"), pt.size=0, ncol=1, cols = c("#cc79a7", "#cc79a7", "#cc79a7", "#cc79a7")
) + NoLegend() + theme(text=element_text(size=30))
p2 | p1
ggsave(args[3], plot = p2 | p1, units = "px", height = 2000, width = 3000)
