# Volcano plot figure 4

library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)

# Load Macosko dataset
macosko.obj <- readRDS(args[1])
macosko.obj <- NormalizeData(macosko.obj)
Idents(macosko.obj) <- "celltype"
macosko.obj$celltype_region <- paste0(macosko.obj$celltype, "_", macosko.obj$regions)
Idents(macosko.obj) <- as.factor(macosko.obj$celltype_region)
idents <- c("astro_CTX", "astro_TH", "neuron_CTX", "neuron_TH")
macosko.obj.ss <- macosko.obj[,macosko.obj$celltype %in% c("neuron", "astro", "Endo", "oligo", "ependymal")]
Idents(macosko.obj.ss) <- as.factor(macosko.obj.ss$celltype)

# Vizgen (default segmentations)
vizgen <- readRDS(args[2])
vizgen.annotated <- readRDS(args[3])
vizgen$annotations <- vizgen.annotated$predicted.class
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
astrocyte.thalamus.cortex.markers <- FindMarkers(vizgen, ident.1 = "Astrocytes_thalamus", ident.2 = "Astrocytes_cortex", logfc.threshold = 0, min.pct = 0, max.cells.per.ident=1000)
astrocyte.thalamus.cortex.markers$gene <- rownames(astrocyte.thalamus.cortex.markers)
genes <- head(rownames(astrocyte.thalamus.cortex.markers), 15)
p <- DotPlot(
  macosko.obj.ss,
  features = genes,
  dot.scale=15,
  split.by = "regions",
  cols = c("blue", "red"))
p <- p + theme(
  text = element_text(size = 45),
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 32),
  axis.text.y = element_text(size = 32))
ggsave(args[8], plot = p, units = "px", height = 4400, width = 6500)

# Vizgen (baysor segmentations)
vizgen.baysor.cort <- readRDS(args[4])
vizgen.baysor.thal <- readRDS(args[5])
vizgen.baysor.thal$celltype_region <- paste0(vizgen.baysor.thal$predicted.class, "_", "thalamus")
vizgen.baysor.cort$celltype_region <- paste0(vizgen.baysor.cort$predicted.class, "_", "cortex")
obj <- merge(vizgen.baysor.cort, vizgen.baysor.thal)
Idents(obj) <- obj$celltype_region
astrocyte.thalamus.cortex.markers <- FindMarkers(obj, ident.1 = "Astrocytes_thalamus", ident.2 = "Astrocytes_cortex", logfc.threshold = 0, min.pct = 0)
astrocyte.thalamus.cortex.markers$gene <- rownames(astrocyte.thalamus.cortex.markers)
genes <- head(rownames(astrocyte.thalamus.cortex.markers), 15)
p <- DotPlot(
  macosko.obj.ss,
  features = genes,
  dot.scale=15,
  split.by = "regions",
  cols = c("blue", "red"))
p <- p + theme(
  text = element_text(size = 45),
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 32),
  axis.text.y = element_text(size = 32))
ggsave(args[9], plot = p, units = "px", height = 4400, width = 6500)

# 10x (baysor segmentations)
tenx.baysor.cort <- readRDS(args[6])
tenx.baysor.thal <- readRDS(args[7])
tenx.baysor.thal$celltype_region <- paste0(tenx.baysor.thal$predicted.class, "_", "thalamus")
tenx.baysor.cort$celltype_region <- paste0(tenx.baysor.cort$predicted.class, "_", "cortex")
obj <- merge(tenx.baysor.cort, tenx.baysor.thal)
Idents(obj) <- obj$celltype_region
astrocyte.thalamus.cortex.markers <- FindMarkers(obj, ident.1 = "Astrocytes_thalamus", ident.2 = "Astrocytes_cortex", logfc.threshold = 0, min.pct = 0)
astrocyte.thalamus.cortex.markers$gene <- rownames(astrocyte.thalamus.cortex.markers)
genes <- head(rownames(astrocyte.thalamus.cortex.markers), 15)
p <- DotPlot(
  macosko.obj.ss,
  features = genes,
  dot.scale=15,
  split.by = "regions",
  cols = c("blue", "red"))
p <- p + theme(
  text = element_text(size = 45),
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 32),
  axis.text.y = element_text(size = 32))
ggsave(args[10], plot = p, units = "px", height = 4400, width = 6500)
