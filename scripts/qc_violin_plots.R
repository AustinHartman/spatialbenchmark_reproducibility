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

# counts per cell
df <- data.frame(
    technology = c(
        rep("Vizgen", ncol(vizgen)),
        rep("10x", ncol(tenx)),
        rep("Resolve", ncol(resolve)),
        rep("EELFISH", ncol(eelfish)),
        rep("MERFISH", ncol(merfish)),
        rep("STARmap PLUS", ncol(starmapplus))
    ),
    counts <- c(
        unlist(vizgen$nCount_RNA),
        unlist(tenx$nCount_RNA),
        unlist(resolve$nCount_RNA),
        unlist(eelfish$nCount_RNA),
        unlist(merfish$nCount_RNA),
        unlist(starmapplus$nCount_RNA)
    )
)
plot <- ggplot(df, aes(x = technology, y = counts, fill = technology)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.25) +
    labs(x = "", y = "") +
    ggtitle("Counts per cell") +
    ylim(0, 750) +
    theme_bw() + scale_fill_manual(values = colors.mapping) + NoLegend() +
    theme(
        axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(args[7], plot, units = "px",  width = 3000, height = 2800)

# features per cell
df <- data.frame(
    technology = c(
        rep("Vizgen", ncol(vizgen)),
        rep("10x", ncol(tenx)),
        rep("Resolve", ncol(resolve)),
        rep("EELFISH", ncol(eelfish)),
        rep("MERFISH", ncol(merfish)),
        rep("STARmap PLUS", ncol(starmapplus))
    ),
    counts <- c(
        unlist(vizgen$nFeature_RNA),
        unlist(tenx$nFeature_RNA),
        unlist(resolve$nFeature_RNA),
        unlist(eelfish$nFeature_RNA),
        unlist(merfish$nFeature_RNA),
        unlist(starmapplus$nFeature_RNA)
    )
)
plot <- ggplot(df, aes(x = technology, y = counts, fill = technology)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.25) +
    labs(x = "", y = "") +
    ggtitle("Features per cell") +
    ylim(0, 200) +
    theme_bw() + scale_fill_manual(values = colors.mapping) + NoLegend() +
    theme(
        axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(args[8], plot, units = "px",  width = 3000, height = 2800)

# counts per feature per cell
df <- data.frame(
    technology = c(
        rep("Vizgen", nrow(vizgen)),
        rep("10x", nrow(tenx)),
        rep("Resolve", nrow(resolve)),
        rep("EELFISH", nrow(eelfish)),
        rep("MERFISH", nrow(merfish)),
        rep("STARmap PLUS", nrow(starmapplus))
    ),
    counts <- c(
        rowSums(vizgen@assays$RNA@counts) / ncol(vizgen),
        rowSums(tenx@assays$RNA@counts) / ncol(tenx),
        rowSums(resolve@assays$RNA@counts) / ncol(resolve),
        rowSums(eelfish@assays$RNA@counts) / ncol(eelfish),
        rowSums(merfish@assays$RNA@counts) / ncol(merfish),
        rowSums(starmapplus@assays$RNA@counts) / ncol(starmapplus))
)
plot <- ggplot(df, aes(x = technology, y = counts, fill = technology)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.25) +
    labs(x = "", y = "") +
    ggtitle("log10(avg. feature counts)") +
    theme_bw() + scale_fill_manual(values = colors.mapping) + NoLegend() + scale_y_log10() +
    theme(
        axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(args[9], plot, units = "px", width = 3000, height = 2800)
