# Generate UMAPs and color cells by different annotations.

suppressWarnings(library(Seurat))
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

vizgen <- readRDS(args[1])
tenx <- readRDS(args[2])
resolve <- readRDS(args[3])
eelfish <- readRDS(args[4])
merfish <- readRDS(args[5])
starmapplus <- readRDS(args[6])
objs <- list(
  vizgen=vizgen, tenx=tenx, resolve=resolve,
  eelfish=eelfish, merfish=merfish, starmapplus=starmapplus)
for (i in names(objs)) {
    objs[[i]] <- NormalizeData(objs[[i]])
    objs[[i]] <- ScaleData(objs[[i]], features = rownames(objs[[i]]))
    objs[[i]] <- RunPCA(objs[[i]], features = rownames(objs[[i]]))
    objs[[i]] <- RunUMAP(objs[[i]], dims = 1:30)
}

# Class (level 1) annotations
custom_colors <- c(
    "Neurons" = "#D32F2F",
    "Neurons,Cycling" = "#1976D2",
    "Astrocyte" = "#388E3C",
    "Ependymal" = "#FBC02D",
    "Immune" = "#8E24AA",
    "OEC" = "#F57C00",
    "Oligos" = "#0288D1",
    "Oligos,Cycling" = "#7B1FA2",
    "Ttr" = "#C2185B",
    "Vascular" = "#7D6608",
    "Satellite-glia" = "#5D4037",
    "Enteric-glia" = "#455A64")
vizgen.plot <- DimPlot(
  objs[["vizgen"]], reduction = "umap", group.by = "predicted.class", raster = FALSE
) + NoLegend() + ggtitle("") + scale_color_manual(values = custom_colors) +
  theme(
    plot.title = element_text(size=25),
    axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)
  ) + xlab("UMAP 1") + ylab("UMAP 2")
tenx.plot <- DimPlot(
  objs[["tenx"]], reduction = "umap", group.by = "predicted.class", raster = FALSE
) + NoLegend() + ggtitle("") + scale_color_manual(values = custom_colors) +
  theme(
    plot.title = element_text(size=25),
    axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)
  ) + xlab("UMAP 1") + ylab("UMAP 2")
resolve.plot <- DimPlot(
  objs[["resolve"]], reduction = "umap", group.by = "predicted.class", raster = FALSE
) + NoLegend() + ggtitle("") + scale_color_manual(values = custom_colors) +
  theme(
    plot.title = element_text(size=25),
    axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)
  ) + xlab("UMAP 1") + ylab("UMAP 2")
eelfish.plot <- DimPlot(
  objs[["eelfish"]], reduction = "umap", group.by = "predicted.class", raster = FALSE
) + NoLegend() + ggtitle("") + scale_color_manual(values = custom_colors) +
  theme(
    plot.title = element_text(size=25),
    axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)
  ) + xlab("UMAP 1") + ylab("UMAP 2")
merfish.plot <- DimPlot(
  objs[["merfish"]], reduction = "umap", group.by = "predicted.class", raster = FALSE
) + NoLegend() + ggtitle("") + scale_color_manual(values = custom_colors) +
  theme(
    plot.title = element_text(size=25),
    axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)
  ) + xlab("UMAP 1") + ylab("UMAP 2")
starmapplus.plot <- DimPlot(
  objs[["starmapplus"]], reduction = "umap", group.by = "predicted.class", raster = FALSE
) + NoLegend() + ggtitle("") + scale_color_manual(values = custom_colors) +
  theme(
    plot.title = element_text(size=25),
    axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)
  ) + xlab("UMAP 1") + ylab("UMAP 2")
plot <- cowplot::plot_grid(vizgen.plot, tenx.plot, resolve.plot, eelfish.plot, merfish.plot, starmapplus.plot, ncol = 3)
ggsave(args[7], plot = plot, width = 24, height = 16)

# Taxonomy group (level 2) annotations
cell_colors <- list(
  "Enteric neurons" = "#4B0082", # Indigo
  "Enteric glia" = "#FF7F50", # Coral
  "Oligodendrocytes" = "#008000", # Green
  "Cholinergic and monoaminergic neurons" = "#9B111E", # Ruby red
  "Telencephalon projecting excitatory neurons" = "#FFA812", # Saffron yellow
  "Telencephalon inhibitory interneurons" = "#B57EDC", # Lavender purple
  "Olfactory inhibitory neurons" = "#40E0D0", # Turquoise
  "Peptidergic neurons" = "#A52A2A", # Brown
  "Di- and mesencephalon excitatory neurons" = "#0000FF", # Blue
  "Glutamatergic neuroblasts" = "#FFD700", # Gold
  "Hindbrain neurons" = "#800080", # Purple
  "Spinal cord excitatory neurons" = "#FF6347", # Tomato
  "Telencephalon projecting inhibitory neurons" = "#20B2AA", # Light Sea Green
  "Olfactory ensheathing cells" = "#FFC0CB", # Pink
  "Non-glutamatergic neuroblasts" = "#3CB371", # Medium Sea Green
  "Dentate gyrus radial glia-like cells" = "#4682B4", # Steel Blue
  "Subventricular zone radial glia-like cells" = "#D2B48C", # Tan
  "Oligodendrocyte precursor cells" = "#8A2BE2", # Blue Violet
  "Ependymal cells" = "#A52A2A", # Brown
  "Subcommissural organ hypendymal cells" = "#DEB887", # Burly Wood
  "Dentate gyrus granule neurons" = "#5F9EA0", # Cadet Blue
  "Cerebellum neurons" = "#D2691E", # Chocolate
  "Di- and mesencephalon inhibitory neurons" = "#FF4500", # Orange Red
  "Spinal cord inhibitory neurons" = "#2E8B57", # Sea Green
  "Vascular and leptomeningeal cells" = "#FF1493", # Deep Pink
  "Vascular smooth muscle cells" = "#1E90FF", # Dodger Blue
  "Pericytes" = "#BDB76B", # Dark Khaki
  "Vascular endothelial cells" = "#8B0000", # Dark Red
  "Microglia" = "#556B2F", # Dark Olive Green
  "Perivascular macrophages" = "#FF00FF", # Magenta
  "Schwann cells" = "#00CED1", # Dark Turquoise
  "Satellite glia" = "#9400D3", # Dark Violet
  "Astrocytes" = "#FFDAB9", # Peach Puff
  "Choroid epithelial cells" = "#ADFF2F", # Green Yellow
  "Peripheral sensory peptidergic neurons" = "#F0E68C", # Khaki
  "Peripheral sensory neurofilament neurons" = "#E9967A", # Dark Salmon
  "Peripheral sensory non-peptidergic neurons" = "#8FBC8F" # Dark Sea Green
)
vizgen.plot <- DimPlot(
  objs[["vizgen"]], reduction = "umap", group.by = "predicted.taxonomy_group", raster = FALSE
) + NoLegend() + ggtitle("") + scale_color_manual(values = cell_colors) +
  theme(
    plot.title = element_text(size=25),
    axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)
  ) + xlab("UMAP 1") + ylab("UMAP 2")
tenx.plot <- DimPlot(
  objs[["tenx"]], reduction = "umap", group.by = "predicted.taxonomy_group", raster = FALSE
) + NoLegend() + ggtitle("") + scale_color_manual(values = cell_colors) +
  theme(
    plot.title = element_text(size=25),
    axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)
  ) + xlab("UMAP 1") + ylab("UMAP 2")
resolve.plot <- DimPlot(
  objs[["resolve"]], reduction = "umap", group.by = "predicted.taxonomy_group", raster = FALSE
) + NoLegend() + ggtitle("") + scale_color_manual(values = cell_colors) +
  theme(
    plot.title = element_text(size=25),
    axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)
  ) + xlab("UMAP 1") + ylab("UMAP 2")
eelfish.plot <- DimPlot(
  objs[["eelfish"]], reduction = "umap", group.by = "predicted.taxonomy_group", raster = FALSE
) + NoLegend() + ggtitle("") + scale_color_manual(values = cell_colors) +
  theme(
    plot.title = element_text(size=25),
    axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)
  ) + xlab("UMAP 1") + ylab("UMAP 2")
merfish.plot <- DimPlot(
  objs[["merfish"]], reduction = "umap", group.by = "predicted.taxonomy_group", raster = FALSE
) + NoLegend() + ggtitle("") + scale_color_manual(values = cell_colors) +
  theme(
    plot.title = element_text(size=25),
    axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)
  ) + xlab("UMAP 1") + ylab("UMAP 2")
starmapplus.plot <- DimPlot(
  objs[["starmapplus"]], reduction = "umap", group.by = "predicted.taxonomy_group", raster = FALSE
) + NoLegend() + ggtitle("") + scale_color_manual(values = cell_colors) +
  theme(
    plot.title = element_text(size=25),
    axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)
  ) + xlab("UMAP 1") + ylab("UMAP 2")
plot <- cowplot::plot_grid(vizgen.plot, tenx.plot, resolve.plot, eelfish.plot, merfish.plot, starmapplus.plot, ncol = 3)
ggsave(args[8], plot = plot, width = 24, height = 16)
