# The purpose of this plot is to show the panel overlap in mouse brain datasets
# across various technologies
# EEL-FISH, 10x, Vizgen, Resolve, osmFISH, starMAP, MERFISH, Resolve

library(Seurat)
library(UpSetR)
library(ggplotify)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

# Vizgen genes
vizgen.genes <- rownames(readRDS(args[1]))

# 10x genes
tenx.genes <- rownames(readRDS(args[2]))

# Resolve
resolve.genes <- rownames(readRDS(args[3]))

# EEL-FISH
eelfish.genes <- rownames(readRDS(args[4]))

# MERFISH 2023 genes
merfish.genes <- rownames(readRDS(args[5]))

# starMAP genes
starmap.genes <- rownames(readRDS(args[6]))

# Gfap is the one shared gene!
all.intersect <- length(
  intersect(eelfish.genes,
    intersect(resolve.genes,
      intersect(vizgen.genes,
        intersect(tenx.genes,
          intersect(starmap.genes, merfish.genes))))))

print(paste("All shared genes:",   intersect(eelfish.genes,
    intersect(resolve.genes,
      intersect(vizgen.genes,
        intersect(tenx.genes,
          intersect(starmap.genes, merfish.genes)))))))

# Dataset
input <- c(
  "10x" = length(tenx.genes),
  "EEL-FISH" = length(eelfish.genes),
  "Vizgen" = length(vizgen.genes),
  "STARmap" = length(starmap.genes),
  "MERFISH" = length(merfish.genes),
  "Resolve" = length(resolve.genes),
  "10x&Vizgen" = length(intersect(tenx.genes, vizgen.genes)),
  "10x&Resolve" = length(intersect(tenx.genes, resolve.genes)),
  "Resolve&Vizgen" = length(intersect(vizgen.genes, resolve.genes)),
  "EEL-FISH&MERFISH" = length(intersect(eelfish.genes, merfish.genes)),
  "STARmap&MERFISH" = length(intersect(starmap.genes, merfish.genes)),
  "EEL-FISH&STARmap" = length(intersect(eelfish.genes, starmap.genes)),
  "10x&Vizgen&Resolve&EEL-FISH&MERFISH&STARmap" = all.intersect
)

# Plot
p <- upset(fromExpression(input),
      nintersects = 40, 
      nsets = 8, 
      order.by = "freq",
      show.numbers = "yes",
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0,
      text.scale = 2.2,
      point.size = 7,
      line.size = 2
)

ggp <- as.ggplot(p, order.by="freq")
ggsave(
  args[7],
  plot = ggp,
  units = "px",
  width = 2600,
  height = 2500)