# Purpose of this script is to generate UMAP plots colored by annotation
# as well as a violin plot showing number of molecules per cell by tech

library(Seurat)
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
print(args)

vizgen <- readRDS(args[1])
tenx <- readRDS(args[3])
resolve <- readRDS(args[5])
eelfish <- readRDS(args[7])
merfish <- readRDS(args[9])
starmapplus <- readRDS(args[11])
objs <- list(
  vizgen=vizgen, tenx=tenx, resolve=resolve,
  eelfish=eelfish, merfish=merfish, starmapplus=starmapplus)

objs[["vizgen"]][["fov"]] <- readRDS(args[2])[["fov"]]
objs[["tenx"]][["fov"]] <- readRDS(args[4])[["fov"]]
objs[["resolve"]][["fov"]] <- readRDS(args[6])[["fov"]]
objs[["eelfish"]][["fov"]] <- readRDS(args[8])[["fov"]]
objs[["merfish"]][["fov"]] <- readRDS(args[10])[["fov"]]
objs[["starmapplus"]][["fov"]] <- readRDS(args[12])[["fov"]]

DefaultBoundary(objs[["vizgen"]][["fov"]]) <- "centroids"
DefaultBoundary(objs[["tenx"]][["fov"]]) <- "centroids"
DefaultBoundary(objs[["resolve"]][["fov"]]) <- "centroids"
DefaultBoundary(objs[["eelfish"]][["fov"]]) <- "centroids"
DefaultBoundary(objs[["merfish"]][["fov"]]) <- "centroids"
DefaultBoundary(objs[["starmapplus"]][["fov"]]) <- "centroids"

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

vizgen.plot <- ImageDimPlot(
  objs[["vizgen"]], group.by = "predicted.class", size = 1,
) + NoLegend() + ggtitle("") + scale_fill_manual(values = custom_colors)
tenx.plot <- ImageDimPlot(
  objs[["tenx"]], group.by = "predicted.class"
) + NoLegend() + ggtitle("") + scale_fill_manual(values = custom_colors)
resolve.plot <- ImageDimPlot(
  objs[["resolve"]], group.by = "predicted.class", size=1,
) + NoLegend() + ggtitle("") + scale_fill_manual(values = custom_colors)
eelfish.plot <- ImageDimPlot(
  objs[["eelfish"]], group.by = "predicted.class"
) + NoLegend() + ggtitle("") + scale_fill_manual(values = custom_colors)
merfish.plot <- ImageDimPlot(
  objs[["merfish"]], group.by = "predicted.class"
) + NoLegend() + ggtitle("") + scale_fill_manual(values = custom_colors)
starmapplus.plot <- ImageDimPlot(
  objs[["starmapplus"]], group.by = "predicted.class"
) + NoLegend() + ggtitle("") + scale_fill_manual(values = custom_colors)

plot <- cowplot::plot_grid(vizgen.plot, tenx.plot, resolve.plot, eelfish.plot, merfish.plot, starmapplus.plot, ncol = 3)
ggsave(args[13], plot = plot, width = 24, height = 16)
