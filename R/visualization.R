library(ggplot2)
library(Seurat)
library(SeuratObject)
library(scales)
library(patchwork)
library(ggrepel)

plot.average.expression <- function(
    obj1,
    obj2,
    group.by.obj1,
    group.by.obj2,
    celltype,
    markers.obj1,
    markers.obj2,
    obj1.name,
    obj2.name,
    non.markers.obj1 = NULL,
    non.markers.obj2 = NULL
) {
    DefaultAssay(obj1) <- "RNA"
    DefaultAssay(obj2) <- "RNA"

    Idents(obj1) <- group.by.obj1
    Idents(obj2) <- group.by.obj2

    obj1.cells <- WhichCells(obj1, idents = celltype)
    obj2.cells <- WhichCells(obj2, idents = celltype)

    overlap.genes <- intersect(rownames(obj1), rownames(obj2))
    message(paste("Some of the shared genes:", head(overlap.genes)))

    obj1.counts <- obj1[["RNA"]]@counts[overlap.genes, obj1.cells]
    obj1.avg <- rowSums(obj1.counts) / ncol(obj1.counts)
    obj2.counts <- obj2[["RNA"]]@counts[overlap.genes, obj2.cells]
    obj2.avg <- rowSums(obj2.counts) / ncol(obj2.counts)
    df <- data.frame("object1" = obj1.avg, "object2" = obj2.avg)
    df$gene <- rownames(df)
    df$marker <- "Non-marker"
    markers.obj1 <- intersect(markers.obj1, overlap.genes)
    markers.obj2 <- intersect(markers.obj2, overlap.genes)
    if (nrow(df[df$gene %in% markers.obj1, ]) > 0) {
        df[df$gene %in% markers.obj1, ]$marker <- obj1.name
    }
    if (nrow(df[df$gene %in% markers.obj2,]) > 0) {
        df[df$gene %in% markers.obj2, ]$marker <- obj2.name
    }
    if (nrow(df[df$gene %in% intersect(markers.obj1, markers.obj2), ]) > 0) {
        df[df$gene %in% intersect(markers.obj1, markers.obj2), ]$marker <- "Marker in both"
    }
    if (!is.null(non.markers.obj1) && !is.null(non.markers.obj2)) {
        df[df$gene %in% union(non.markers.obj1, non.markers.obj2), ]$marker <- "Other Celltype marker"
    }
    p1 <- ggplot(
    df, mapping = aes(x = object2, y = object1, color = marker)
    ) + geom_point(size = 2.2) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    geom_hline(yintercept = 1) +
    geom_vline(xintercept = 1) +
    scale_x_continuous(
      trans = "log2",
      breaks = trans_breaks("log2", function(x) 2^x),
      labels = trans_format("log2", math_format(2^.x))
    ) +
    scale_y_continuous(
      trans = "log2",
      breaks = trans_breaks("log2", function(x) 2^x),
      labels = trans_format("log2", math_format(2^.x))
    ) +
    theme_bw() + annotation_logticks() +
    geom_text_repel(aes(label = ifelse(marker == "bothmarker", gene, "")), color = "black") +
    ylab(obj1.name) + xlab(obj2.name) +
    theme(text = element_text(size = 25))
    title <- paste0("Average ", celltype, " Gene Expression")
    color.map <- c(
        "bothmarker" = "red",
        "nonmarker" = "grey85",
        "obj1marker" = "blue",
        "obj2marker" = "darkgreen",
        "Other Cell Marker" = "grey65")
    names(color.map) <- c("Marker in both", "Non-marker", obj1.name, obj2.name, "Other Celltype marker")
    return(p1 + ggtitle(title) + scale_color_manual(values = color.map))
}

plot.feature.scatters <- function(
    obj,
    gene1,
    gene2,
    color,
    max.x = 100,
    max.y = 100,
    size = 0.2
) {
    data <- FetchData(obj, vars = c(gene1, gene2), slot = "counts")
    colnames(data) <- c("gene1", "gene2")
    title <- paste(gene1, "vs.", gene2)
    p <- ggplot(data=data, mapping=aes(x=gene1, y=gene2)) +
        geom_point(alpha=0.5, size=size, color = color) +
        xlim(0, max.x) + ylim(0, max.y) + ggtitle(title) + theme_classic() +
        theme(
        plot.title = element_text(size=25),
        axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)
        ) + xlab(gene1) + ylab(gene2)

    return(p)
}

marker.heatmap <- function(
    object,
    markers,
    group.by
) {
  Idents(object) <- group.by
  heatmap.genes <- c()
  for (ct in unique(markers$cluster)) {
    heatmap.genes <- c(
      heatmap.genes,
      head(markers[markers$cluster == ct,]$gene, 10)
    )
  }
  cells <- c()
  for (ct in unique(markers$cluster)) {
    ct.cells <- WhichCells(object, idents = ct)
    if (length(ct.cells) > 400) {
      cells <- c(cells, sample(ct.cells, 400))
    } else {
      cells <- c(cells, ct.cells)
    }
  }
  object <- ScaleData(object, features = heatmap.genes)
  object.heatmap <- DoHeatmap(
    object, features = heatmap.genes, cells = cells, disp.max = 2
  ) + ggtitle("Heatmap")
  return(object.heatmap)
}
