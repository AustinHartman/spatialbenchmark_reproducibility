# Volcano plots supplementary figure 6

library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)

### Xenium ###
tenx.baysor.cort <- readRDS(args[1])
tenx.baysor.thal <- readRDS(args[2])
tenx.baysor.thal$celltype_region <- paste0(tenx.baysor.thal$predicted.class, "_", "thalamus")
tenx.baysor.cort$celltype_region <- paste0(tenx.baysor.cort$predicted.class, "_", "cortex")
obj <- merge(tenx.baysor.cort, tenx.baysor.thal)
Idents(obj) <- obj$celltype_region
astrocyte.thalamus.cortex.markers <- FindMarkers(
  obj, ident.1 = "Astrocytes_thalamus", ident.2 = "Astrocytes_cortex",
  logfc.threshold = 0, min.pct = 0)
astrocyte.thalamus.cortex.markers$gene <- rownames(astrocyte.thalamus.cortex.markers)
p <- ggplot(astrocyte.thalamus.cortex.markers, aes(
  avg_log2FC,
  -log10(p_val),
  label=ifelse(gene %in% head(rownames(astrocyte.thalamus.cortex.markers), 15), gene, "")
)
) +
  geom_point(size = 3) + geom_text_repel(size = 8) + theme_minimal() + theme(
    text = element_text(size=30),
    axis.line.x=element_line(size=2),
    axis.line.y=element_line(size=2)
  ) + xlab("") + ylab("")
ggsave(args[5], plot = p, units = "px", height = 2000, width = 2400)
Idents(obj) <- as.factor(obj$celltype_region)
obj <- NormalizeData(obj)
p1<-VlnPlot(
  obj, idents=c("Astrocytes_cortex", "Astrocytes_thalamus", "Neurons_cortex", "Neurons_thalamus"),
  cols=c("#d55e00","#d55e00","#d55e00","#d55e00"), features=c("Slc17a6"), pt.size=0
) + theme(text=element_text(size = 34)) + ylab("") + xlab("") + NoLegend()
p2<-VlnPlot(
  obj, idents=c("Astrocytes_cortex", "Astrocytes_thalamus", "Neurons_cortex", "Neurons_thalamus"),
  cols=c("#cc79a7","#cc79a7","#cc79a7","#cc79a7"), features=c("Satb2"), pt.size=0
) + theme(text=element_text(size = 34)) + ylab("") + xlab("") + NoLegend()
ggsave(args[6], plot = p2 | p1, units = "px", height = 2000, width = 3000)

### MERSCOPE ###
vizgen.baysor.cort <- readRDS(args[3])
vizgen.baysor.thal <- readRDS(args[4])
vizgen.baysor.thal$celltype_region <- paste0(vizgen.baysor.thal$predicted.class, "_", "thalamus")
vizgen.baysor.cort$celltype_region <- paste0(vizgen.baysor.cort$predicted.class, "_", "cortex")
obj <- merge(vizgen.baysor.cort, vizgen.baysor.thal)
Idents(obj) <- obj$celltype_region
astrocyte.thalamus.cortex.markers <- FindMarkers(obj, ident.1 = "Astrocytes_thalamus", ident.2 = "Astrocytes_cortex", logfc.threshold = 0, min.pct = 0)
astrocyte.thalamus.cortex.markers$gene <- rownames(astrocyte.thalamus.cortex.markers)
p <- ggplot(astrocyte.thalamus.cortex.markers, aes(
  avg_log2FC,
  -log10(p_val),
  label=ifelse(gene %in% head(rownames(astrocyte.thalamus.cortex.markers), 15), gene, "")
)
) +
  geom_point(size = 3) + geom_text_repel(size = 8) + theme_minimal() + theme(
    text = element_text(size=30),
    axis.line.x=element_line(size=2),
    axis.line.y=element_line(size=2)
  ) + xlab("") + ylab("")
ggsave(args[7], plot = p, units = "px", height = 2000, width = 2400)
Idents(obj) <- as.factor(obj$celltype_region)
obj <- NormalizeData(obj)
p1<-VlnPlot(
  obj, idents=c("Astrocytes_cortex", "Astrocytes_thalamus", "Neurons_cortex", "Neurons_thalamus"),
  cols=c("#d55e00","#d55e00","#d55e00","#d55e00"), features=c("Slc17a6"), pt.size=0
) + theme(text=element_text(size = 34)) + ylab("") + xlab("") + NoLegend()
p2<-VlnPlot(
  obj, idents=c("Astrocytes_cortex", "Astrocytes_thalamus", "Neurons_cortex", "Neurons_thalamus"),
  cols=c("#cc79a7","#cc79a7","#cc79a7","#cc79a7"), features=c("Slc17a7"), pt.size=0
) + theme(text=element_text(size = 34)) + ylab("") + xlab("") + NoLegend()
ggsave(args[8], plot = p2 | p1, units = "px", height = 2000, width = 3000)
