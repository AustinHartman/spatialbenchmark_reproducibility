library(Seurat)
library(future)
library(pbapply)
library(RImageJROI)
library(hdf5r)
library(rlang)

plan("multisession", workers = 24)
source("./R/loadData.R")

args <- commandArgs(trailingOnly = TRUE)
print(args)

result <- switch(
    args[2],
    "vizgen" = LoadVizgen(args[1], fov = "fov", assay = "RNA"),
    "tenx" = LoadXenium(args[1], fov = "fov", assay = "RNA"),
    "merfish" = LoadMERFISH23(args[1], s = args[length(args)]),
    "resolve" = LoadResolve(args[1]),
    "eelfish" = LoadEELFISH(args[1]),
    "starmapplus" = LoadSTARmapPlus(args[1]),
)

if (length(args) > 2) {
    save.path <- file.path("./data/", paste0(args[3], ".rds"))
} else {
    save.path <- file.path("./data/", paste0(args[2], ".rds"))
}
print(save.path)
saveRDS(result, save.path)
