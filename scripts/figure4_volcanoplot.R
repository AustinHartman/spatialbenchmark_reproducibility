# Volcano plot figure 4

library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)

tenx <- readRDS(args[1])
tenx.annotated <- readRDS(args[2])
tenx$annotations <- tenx.annotated$predicted.class

colors <- list(
  "Neurons" = "lightgrey",
  "Astrocytes" = "red",
  "Vascular" = "lightgrey",
  "Immune" = "lightgrey",
  "Oligos" = "lightgrey",
  "Ependymal" = "lightgrey")

# Rename the cells of interest in regions of interest
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

tenx$annotations <- ifelse(
    tenx$annotations == "Astrocytes" & Cells(tenx) %in% thalamus.cells$cell,
    "Astrocytes_thalamus", tenx$annotations)
tenx$annotations <- ifelse(
    tenx$annotations == "Astrocytes" & Cells(tenx) %in% cortex.cells$cell,
    "Astrocytes_cortex", tenx$annotations)

Idents(tenx) <- tenx$annotations
astro.thal.ctx.markers <- FindMarkers(
    tenx, ident.1 = "Astrocytes_thalamus", ident.2 = "Astrocytes_cortex",
    logfc.threshold = 0, min.pct = 0)
print(head(astro.thal.ctx.markers, 70))
print(astro.thal.ctx.markers["Slc17a7", ])
astro.thal.ctx.markers$gene <- rownames(astro.thal.ctx.markers)
p <- ggplot(astro.thal.ctx.markers, aes(
    avg_log2FC,
    -log10(p_val),
    label=ifelse(
        gene %in% head(rownames(astro.thal.ctx.markers), 20),
        gene, "")
  )
) +
  geom_point(size = 3) + geom_text_repel(size = 8, max.overlaps = 1000) +
  theme_minimal() + theme(
    text = element_text(size=30),
    axis.line.x=element_line(size=2),
    axis.line.y=element_line(size=2)
  ) + xlim(-2.5, 2.5) + xlab("") + ylab("")
ggsave(args[3], plot = p, units = "px", height = 2000, width = 2400)
