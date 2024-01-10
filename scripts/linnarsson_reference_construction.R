# Script to load in the Linnarsson mouse reference from the provided loom file

#  Qeustions to answer: how do we integrate this data for annotation across spatial methods?

library(Seurat)
library(hdf5r)

args <- commandArgs(trailingOnly = TRUE)

path <- args[1]

file.h5 <- H5File$new(path, mode="r")
matrix <- file.h5[["matrix"]][,]
cell.ids <- file.h5[["col_attrs"]][["CellID"]][]
gene.names <- file.h5[["row_attrs"]][["Gene"]][]

# load in the dataset/batch information
donor.id <- file.h5[["col_attrs"]][["DonorID"]][]
sample.id <- file.h5[["col_attrs"]][["SampleID"]][]

# load in the cell type information
cell.class <- file.h5[["col_attrs"]][["Class"]][]
cell.subclass <- file.h5[["col_attrs"]][["Subclass"]][]
cell.taxonomy.group <- file.h5[["col_attrs"]][["Taxonomy_group"]][]
cell.taxonomy.three <- file.h5[["col_attrs"]][["TaxonomyRank3"]][]
cell.description <- file.h5[["col_attrs"]][["Description"]][]
cell.sex <- file.h5[["col_attrs"]][["Sex"]][]

rownames(matrix) <- cell.ids
colnames(matrix) <- gene.names

object <- CreateSeuratObject(counts = t(matrix))
object <- AddMetaData(object, donor.id, col.name="DonorID")
object <- AddMetaData(object, sample.id, col.name="SampleID")
object <- AddMetaData(object, cell.class, col.name="Class")
object <- AddMetaData(object, cell.subclass, col.name="Subclass")
object <- AddMetaData(object, cell.description, col.name="Description")
object <- AddMetaData(object, cell.sex, col.name = "Sex")
object <- AddMetaData(object, cell.taxonomy.group, col.name = "Taxonomy_group")
object <- AddMetaData(object, cell.taxonomy.three, col.name = "TaxonomyRank3")

object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^mt-")
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# filter cells with low UMI counts
object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)

# filter samples with low numbers of cells
keep.samples <- names(table(object$SampleID))[table(object$SampleID) > 500]
cells.keep <- Cells(object)[object$SampleID %in% keep.samples]
object <- subset(object, cells = cells.keep)

# unintegrated analysis of the linnarsson dataset to see how it all looks
object <- NormalizeData(object)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 4000)
object <- ScaleData(object, features = VariableFeatures(object))
object <- RunPCA(object, npcs = 100, features = VariableFeatures(object))
object <- RunUMAP(object, dims = 1:100, return.model = TRUE)
object <- FindNeighbors(object, reduction = "pca", dims = 1:100)

saveRDS(object, args[2])
