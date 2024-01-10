args <- commandArgs(trailingOnly = TRUE)
obj <- readRDS(args[1])
write.csv(obj@meta.data, args[2])
