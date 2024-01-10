# The purpose of this script is to generate summary violin plots using each of the technologies

library(Seurat)
library(ggplot2)
library(Matrix)
source("./R/utils.R")

args <- commandArgs(trailingOnly = TRUE)

vizgen <- readRDS(args[1])
tenx <- readRDS(args[2])
resolve <- readRDS(args[3])
eelfish <- readRDS(args[4])
merfish <- readRDS(args[5])
starmapplus <- readRDS(args[6])


technologies <- c("Vizgen", "MERFISH", "Resolve", "STARmap PLUS", "10x", "EELFISH")
colors.mapping <- colors[tolower(technologies)]
names(colors.mapping) <- technologies
colors.mapping <- lapply(colors.mapping, function(x) { ifelse(is.null(x), "black", x) })
print(colors.mapping)

df <- data.frame(
    technology <- technologies,
    panel_size = c(
        nrow(vizgen),
        nrow(merfish),
        nrow(resolve),
        nrow(starmapplus),
        nrow(tenx),
        nrow(eelfish)
    ),
    mols_per_cell = c(
        mean(vizgen$nCount_RNA),
        mean(merfish$nCount_RNA),
        mean(resolve$nCount_RNA),
        mean(starmapplus$nCount_RNA),
        mean(tenx$nCount_RNA),
        mean(eelfish$nCount_RNA)
    )
)
plot <- ggplot(df, aes(x = panel_size, y = mols_per_cell, color = technology)) +
    geom_point(size = 3) +
    # labs(x = "", y = "") +
    ggtitle("Panel size vs Molecules per cell") +
    theme_bw() + scale_color_manual(values = colors.mapping) + NoLegend() +
    theme(
        axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(args[7], plot, units = "px",  width = 3000, height = 2800)

# counts per feature per cell
df <- data.frame(
    technology <- technologies,
    panel_size = c(
        nrow(vizgen),
        nrow(merfish),
        nrow(resolve),
        nrow(starmapplus),
        nrow(tenx),
        nrow(eelfish)
    ),
    mols_per_gene <- c(
        mean(rowSums(vizgen@assays$RNA@counts) / ncol(vizgen)),
        mean(rowSums(merfish@assays$RNA@counts) / ncol(merfish)),
        mean(rowSums(resolve@assays$RNA@counts) / ncol(resolve)),
        mean(rowSums(starmapplus@assays$RNA@counts) / ncol(starmapplus)),
        mean(rowSums(tenx@assays$RNA@counts) / ncol(tenx)),
        mean(rowSums(eelfish@assays$RNA@counts) / ncol(eelfish))
    )
)
plot <- ggplot(df, aes(x = panel_size, y = mols_per_gene, color = technology)) +
    geom_point(size = 3) +
    # labs(x = "", y = "") +
    ggtitle("Panel size vs Molecules per gene") +
    theme_bw() + scale_color_manual(values = colors.mapping) + NoLegend() +
    theme(
        axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(args[8], plot, units = "px",  width = 3000, height = 2800)
