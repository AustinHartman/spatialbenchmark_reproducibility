library(Seurat)

annotate.seurat.objs <- function(
    obj,
    reference,
    refdata = list(
      "class" = "Class",
      "subclass" = "Subclass",
      "description" = "Description",
      "taxonomy_group" = "Taxonomy_group",
      "taxonomy_3" = "TaxonomyRank3"
    )
) {
  # drop images from the object
  if ("fov" %in% Images(obj)) {
    obj[["fov"]] <- NULL
  }

  # annotate vizgen object
  anchors <- FindTransferAnchors(
    reference = reference,
    query = obj,
    reduction = "cca",
    query.assay = "RNA",
    reference.assay = "RNA")
  obj <- MapQuery(
    anchorset = anchors,
    reference = reference,
    query = obj,
    refdata = refdata,
    reduction.model = "umap")
  return(obj)
}

args <- commandArgs(trailingOnly = TRUE)
obj <- readRDS(file.path(getwd(), args[2]))
obj <- NormalizeData(obj)
reference <- readRDS(args[1])
result <- annotate.seurat.objs(obj = obj, reference = reference)
saveRDS(result, args[3])
