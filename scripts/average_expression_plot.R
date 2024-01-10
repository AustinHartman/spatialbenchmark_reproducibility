library(Seurat)
library(SeuratWrappers)
library(stringr)
source("./R/visualization.R")

args <- commandArgs(trailingOnly = TRUE)
options(str = list(strict.width = "no"))
message(paste0("Arguments: ", str(args)))
if (file.exists(args[1])) {
    object <- readRDS(args[1])
} else {
    object <- readRDS(file.path(getwd(), args[1]))
}
reference <- readRDS(args[2])
print(head(rownames(object)))
features <- intersect(rownames(object), rownames(reference))

Idents(object) <- "predicted.subclass"
Idents(reference) <- "Subclass"

reference.markers <- RunPrestoAll(reference, slot = "data", features = features)
object.markers <- RunPrestoAll(object, slot = "data", features = features)

celltype <- args[4]

markers.obj1 <- object.markers[
    object.markers$cluster == celltype & object.markers$avg_log2FC > 0.25,
]
markers.obj2 <- reference.markers[
    reference.markers$cluster == celltype & reference.markers$avg_log2FC > 0.25,
]
non.markers.obj1 <- object.markers[
    object.markers$cluster != celltype & object.markers$avg_log2FC > 0.25,
]$gene
non.markers.obj2 <- reference.markers[
    reference.markers$cluster != celltype & reference.markers$avg_log2FC > 0.25,
]$gene

if (is.null(markers.obj1) || nrow(markers.obj1) == 0) {
    markers.obj1 <- c()
} else {
    markers.obj1 <- markers.obj1$gene
}
if (is.null(markers.obj2) || nrow(markers.obj2) == 0) {
    markers.obj2 <- c()
} else {
    markers.obj2 <- markers.obj2$gene
}

plot <- plot.average.expression(
    obj1 = object,
    obj2 = reference,
    group.by.obj1 = "predicted.subclass",
    group.by.obj2 = "Subclass",
    celltype = celltype,
    markers.obj1 = markers.obj1,
    markers.obj2 = markers.obj2,
    obj1.name = str_to_title(args[3]),
    obj2.name = "Linnarsson Reference",
    non.markers.obj1 = setdiff(non.markers.obj1, markers.obj1),
    non.markers.obj2 = setdiff(non.markers.obj2, markers.obj2)
) + NoLegend()

if (length(args) >= 5 && args[5] == "FALSE") {
    plot <- plot + scale_color_manual(
        values = c("black", "black", "black", "black", "black"))
    save.path <- file.path(
    "./data/",
    paste0("average_expression_", args[4], "_", args[3], "_no_color.png"))
    ggplot2::ggsave(
    filename = save.path,
    plot = plot,
    units = "px",
    width = 3000,
    height = 2800)
} else {
    save.path <- file.path(
    "./data/",
    paste0("average_expression_", args[4], "_", args[3], ".png"))
    ggplot2::ggsave(
    filename = save.path,
    plot = plot,
    units = "px",
    width = 2800,
    height = 2800)
}
