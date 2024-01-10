# The purpose of this plot is to generate counts per feature rank plots for
# each of the available technologies with background probes highlighted.

library(Seurat)
library(ggplot2)
library(scales)

args <- commandArgs(trailingOnly = TRUE)
molecules <- read.csv(args[1])
technology <- args[2]

color.mapping <- list(
    "Gene probes" = "black",
    "Negative probes" = "firebrick1")
counts <- table(molecules$gene)

if (technology == "vizgen") {
    background.features <- names(counts)[grepl("Blank*", names(counts))]
} else if (technology == "tenx") {
    background.features <- names(counts)[grepl("NegControl*|BLANK*", names(counts))]
} else if (technology == "resolve") {
    background.features <- names(counts)[grepl("FP*", names(counts))]
} else if (technology == "eelfish") {
    background.features <- names(counts)[grepl("Control*", names(counts))]
} else if (technology == "merfish") {
    background.features <- names(counts)[grepl("blank*", names(counts))]
}

df <- data.frame(
    features = names(counts),
    counts = unlist(counts),
    category = "Gene probes")
df[df$features %in% background.features, "category"] <- "Negative probes"
write.csv(df, paste0("/brahms/hartmana/spatial_sensitivity_comparison/spatial_benchmarking_rank_counts_", technology))
plot <- ggplot(df, aes(x = reorder(features, -counts), y = counts, fill = category)) +
    geom_col(stat = "identity", width = 1) + scale_y_continuous(labels = label_comma()) +
    labs(x = "Feature", y = "Counts") + NoLegend() +
    theme_classic() + scale_fill_manual(values = color.mapping) +
    theme(axis.text = element_text(color = "black", size = 30), axis.title = element_text(size = 28), axis.text.x = element_blank())
ggsave(args[length(args) - 1], plot, units = "px",  width = 3800, height = 3400)
if (technology == "resolve") {
    ggsave(args[length(args)], plot + ylim(0, 50000), units = "px",  width = 3800, height = 3400)
} else {
    ggsave(args[length(args)], plot + ylim(0, 15000), units = "px",  width = 3800, height = 3400)
}
