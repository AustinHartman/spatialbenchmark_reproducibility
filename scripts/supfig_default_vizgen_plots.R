# Volcano plot figure 4

library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
vizgen <- readRDS(args[1])
vizgen.annotated <- readRDS(args[2])
vizgen$annotations <- vizgen.annotated$predicted.class

colors <- list(
  "Neurons" = "lightgrey",
  "Astrocytes" = "red",
  "Vascular" = "lightgrey",
  "Immune" = "lightgrey",
  "Oligos" = "lightgrey",
  "Ependymal" = "lightgrey"
)

cell.names <- vizgen@images$fov@boundaries$centroids@cells
pos <- vizgen@images$fov@boundaries$centroids@coords
df <- as.data.frame(pos)
df$cell <- cell.names

thalamus.cells <- df[
  (df$x>3250 & df$x<4000 & df$y>3250 & df$y<4000)
  , ]
cortex.cells <- df[
  (df$x>3500 & df$x<4500 & df$y>500 & df$y<1500)
  , ]
vizgen$annotations <- ifelse(vizgen$annotations == "Astrocytes" & Cells(vizgen) %in% thalamus.cells$cell, "Astrocytes_thalamus", vizgen$annotations)
vizgen$annotations <- ifelse(vizgen$annotations == "Astrocytes" & Cells(vizgen) %in% cortex.cells$cell, "Astrocytes_cortex", vizgen$annotations)
vizgen$annotations <- ifelse(vizgen$annotations == "Neurons" & Cells(vizgen) %in% thalamus.cells$cell, "Neurons_thalamus", vizgen$annotations)
vizgen$annotations <- ifelse(vizgen$annotations == "Neurons" & Cells(vizgen) %in% cortex.cells$cell, "Neurons_cortex", vizgen$annotations)
Idents(vizgen) <- vizgen$annotations
astrocyte.thalamus.cortex.markers <- FindMarkers(vizgen, ident.1 = "Astrocytes_thalamus", ident.2 = "Astrocytes_cortex", logfc.threshold = 0, min.pct = 0)
astrocyte.thalamus.cortex.markers$gene <- rownames(astrocyte.thalamus.cortex.markers)
p <- ggplot(astrocyte.thalamus.cortex.markers, aes(
  avg_log2FC,
  -log10(p_val),
  label=ifelse(gene %in% head(rownames(astrocyte.thalamus.cortex.markers), 15), gene, "")
)
) +
  geom_point(size = 3) + geom_text_repel(size = 8) + theme_minimal() + theme(
    text = element_text(size=30),
    axis.line.x=element_line(size=2),
    axis.line.y=element_line(size=2)
  ) + xlab("") + ylab("")
ggsave(args[3], plot = p, units = "px", height = 2000, width = 2400)

Idents(vizgen) <- as.factor(vizgen$annotations)
vizgen <- NormalizeData(vizgen)
p1<-VlnPlot(
  vizgen, idents=c("Astrocytes_cortex", "Astrocytes_thalamus", "Neurons_cortex", "Neurons_thalamus"),
  cols=c("#d55e00","#d55e00","#d55e00","#d55e00"), features=c("Slc17a6"), pt.size=0
) + theme(text=element_text(size = 34)) + ylab("") + xlab("") + NoLegend()
p2<-VlnPlot(
  vizgen, idents=c("Astrocytes_cortex", "Astrocytes_thalamus", "Neurons_cortex", "Neurons_thalamus"),
  cols=c("#cc79a7","#cc79a7","#cc79a7","#cc79a7"), features=c("Slc17a7"), pt.size=0
) + theme(text=element_text(size = 34)) + ylab("") + xlab("") + NoLegend()
ggsave(args[4], plot = p2 | p1, units = "px", height = 2000, width = 3000)
