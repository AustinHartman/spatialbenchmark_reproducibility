# The purpose of this plot is to show the panel overlap in mouse brain datasets
# across various technologies
# EEL-FISH, 10x, Vizgen, Resolve, osmFISH, starMAP, MERFISH, Resolve

library(Seurat)
library(UpSetR)
library(ggplotify)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

vizgen <- readRDS(args[1])
tenx <- readRDS(args[2])
resolve <- readRDS(args[3])
eelfish <- readRDS(args[4])
merfish <- readRDS(args[5])
starmap <- readRDS(args[6])

df <- data.frame(
    tech = c("Vizgen", "10x", "Resolve", "EEL-FISH", "MERFISH", "STARmap PLUS"),
    avgcounts = c(mean(vizgen$nCount_RNA), mean(tenx$nCount_RNA), mean(resolve$nCount_RNA), mean(eelfish$nCount_RNA), mean(merfish$nCount_RNA), mean(starmap$nCount_RNA)),
    ngenes = c(nrow(vizgen), nrow(tenx), nrow(resolve), nrow(eelfish), nrow(merfish), nrow(starmap))
)

p <- ggplot(df, aes(y = avgcounts, x = ngenes)) + geom_point(size = 5) + theme_bw() + theme(text = element_text(size = 30, color = "black"))
ggsave(
  args[7],
  plot = p,
  units = "px",
  width = 3200,
  height = 2500)