# Save an exclusive expression plot

library(ggplot2) 
library(Seurat)
source("./R/visualization.R")
source("./R/utils.R")

args <- commandArgs(trailingOnly = TRUE)
message(paste("Arguments: ", args))

max.x <- ifelse(length(args) >= 5, as.numeric(args[5]), 100)
max.y <- ifelse(length(args) >= 6, as.numeric(args[6]), 100)
color <- ifelse(length(args) > 7, args[7], colors[tolower(args[2])])
color <- ifelse(is.null(color), "black", color)
size <- ifelse(length(args) >= 8, as.numeric(args[8]), 0.5)

if (file.exists(args[1])) {
    obj <- readRDS(args[1])
} else {
    obj <- readRDS(file.path(getwd(), args[1]))
}

plot <- plot.feature.scatters(
    obj = obj,
    gene1 = args[3],
    gene2 = args[4],
    color = color,
    max.x = max.x,
    max.y = max.y,
    size = size
)
ggsave(filename = args[length(args)], plot = plot, units = "px", width = 1500, height = 1500)
