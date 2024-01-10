# scRNA-seq regional cell type specific expression plots

library(Seurat)
library(ggplot2)
library(SeuratWrappers)

args <- commandArgs(trailingOnly = TRUE)
macosko.obj <- readRDS(args[1])
macosko.obj <- NormalizeData(macosko.obj)
Idents(macosko.obj) <- "celltype"
macosko.obj$celltype_region <- paste0(macosko.obj$celltype, "_", macosko.obj$regions)
Idents(macosko.obj) <- as.factor(macosko.obj$celltype_region)
idents <- c("astro_CTX", "astro_TH", "neuron_CTX", "neuron_TH")

p1 <- VlnPlot(macosko.obj, features = "Satb2", pt.size = 0, idents=idents, cols=c("#cc79a7","#cc79a7","#cc79a7","#cc79a7")) + NoLegend() + ggtitle("Satb2") + theme(text=element_text(size=30))
p2 <- VlnPlot(macosko.obj, features = "Slc17a6", pt.size = 0, idents=idents, cols=c("#d55e00","#d55e00","#d55e00","#d55e00")) + NoLegend() + ggtitle("Slc17a6") + theme(text=element_text(size=30))
p1 | p2
ggsave(args[2], plot = p1 | p2, units = "px", height = 2000, width = 3000)

macosko.obj.ss <- macosko.obj[,macosko.obj$celltype %in% c("neuron", "astro", "Endo", "oligo", "ependymal")]
Idents(macosko.obj.ss) <- as.factor(macosko.obj.ss$celltype)
markers <- c(
  "Satb2", "Lamp5", "Dkk3", "Neurod6", "Arc", "Fezf2",
  "Slc17a6", "Rnf152", "Unc13c", "Dner", "Prox1", "Nr2f2", "Necab2",
  "Hapln1", "Clmn", "Meis2",
  "Tanc1", "Rims3", "Pde7b")
p <- DotPlot(
  macosko.obj.ss,
  features = markers,
  dot.scale=15,
  split.by = "regions",
  cols = c("blue", "red"))
p <- p + theme(
  text = element_text(size = 45),
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 32),
  axis.text.y = element_text(size = 32))
ggsave(args[3], plot = p, units = "px", height = 4400, width = 6500)
