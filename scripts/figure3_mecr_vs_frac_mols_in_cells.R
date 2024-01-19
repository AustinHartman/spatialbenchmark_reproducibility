# Construct a scatter plot displaying the MECR vs. the fraction of molecules

library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(Matrix)
source("./R/utils.R")

technologies <- c("Vizgen", "MERFISH", "Resolve", "10x", "EELFISH")
colors.mapping <- colors[tolower(technologies)]
names(colors.mapping) <- technologies
colors.mapping <- lapply(colors.mapping, function(x) { ifelse(is.null(x), "black", x) })

args <- commandArgs(trailingOnly = TRUE)

######################################
### Fraction of molecules in cells ###
######################################

vizgen <- readRDS(args[1])
tenx <- readRDS(args[2])
resolve <- readRDS(args[3])
eelfish <- readRDS(args[4])
merfish <- readRDS(args[5])
message("Seurat objects loaded")

vizgen.mols <- read.csv(args[6])
vizgen.genes <- vizgen.mols$gene
tenx.mols <- read.csv(args[7])
tenx.genes <- tenx.mols$gene
resolve.mols <- read.csv(args[8])
resolve.genes <- resolve.mols$gene
eelfish.mols <- read.csv(args[9])
eelfish.genes <- eelfish.mols$gene
merfish.mols <- read.csv(args[10])
merfish.genes <- merfish.mols$gene
message("Molecules loaded")

# filter names of negative probes
vizgen.bgd.feats <- names(vizgen.genes)[grepl("Blank*", names(vizgen.genes))]
tenx.bgd.feats <- names(tenx.genes)[grepl("NegControl*|BLANK*", names(tenx.genes))]
resolve.bgd.feats <- names(resolve.genes)[grepl("^FP.*", names(resolve.genes))]
eelfish.bgd.feats <- names(eelfish.genes)[grepl("Control*", names(eelfish.genes))]
merfish.bgd.feats <- names(merfish.genes)[grepl("blank*", names(merfish.genes))]
message("Negative probes identified")

# remove negative probes from counts
vizgen.counts <- vizgen@assays$RNA@counts[
    !(rownames(vizgen@assays$RNA@counts) %in% vizgen.bgd.feats), ]
tenx.counts <- tenx@assays$RNA@counts[
    !(rownames(tenx@assays$RNA@counts) %in% tenx.bgd.feats), ]
resolve.counts <- resolve@assays$RNA@counts[
    !(rownames(resolve@assays$RNA@counts) %in% resolve.bgd.feats), ]
eelfish.counts <- eelfish@assays$RNA@counts[
    !(rownames(eelfish@assays$RNA@counts) %in% eelfish.bgd.feats), ]
merfish.counts <- merfish@assays$RNA@counts[
    !(rownames(merfish@assays$RNA@counts) %in% merfish.bgd.feats), ]
message("Negative probes filtered from counts")

# remove negative probe molecules
vizgen.mols <- vizgen.mols[!(vizgen.mols$gene %in% vizgen.bgd.feats), ]
tenx.mols <- tenx.mols[!(tenx.mols$gene %in% tenx.bgd.feats), ]
resolve.mols <- resolve.mols[!(resolve.mols$gene %in% resolve.bgd.feats), ]
eelfish.mols <- eelfish.mols[!(eelfish.mols$gene %in% eelfish.bgd.feats), ]
merfish.mols <- merfish.mols[!(merfish.mols$gene %in% merfish.bgd.feats), ]
message("Negative probes filtered from molecules")

# total number of molecules in each dataset
# (including molecules unassigned to a cell)
vizgen.nmols <- nrow(vizgen.mols)
tenx.nmols <- nrow(tenx.mols)
resolve.nmols <- nrow(resolve.mols)
eelfish.nmols <- nrow(eelfish.mols)
merfish.nmols <- nrow(merfish.mols)

########################
### MECR calculation ###
########################

# load in other objects for comparison (not-necessarily spatial)
objects <- list(
    "Vizgen" = readRDS(args[1]),
    "10x" = readRDS(args[2]),
    "Resolve" = readRDS(args[3]),
    "MERFISH" = readRDS(args[4]),
    "EELFISH" = readRDS(args[5]))

# identify markers in the reference (Linnarsson)
ref <- readRDS(args[length(args) - 1])
ref <- NormalizeData(ref, verbose = FALSE)
Idents(ref) <- "Class"
ref.markers <- RunPrestoAll(ref)
ref.markers.ss <- ref.markers[
    ref.markers$pct.1 > 0.25 & ref.markers$pct.2 < 0.01, ]

get.coexpression.rate <- function(obj, sc.markers) {
    coexp.rates <- c()
    genes <- intersect(rownames(obj), rownames(sc.markers))
    print(length(genes))
    if (length(genes) > 25) genes <- sample(genes, 25)
    mtx <- as.matrix(obj@assays$RNA@counts[genes, ])
    for (g1 in genes) {
        for (g2 in genes) {
            if (
                (g1 != g2) && (g1 > g2) &&
                (sc.markers[g1, "cluster"] != sc.markers[g2, "cluster"])
            ) {
                c1 <- mtx[g1, ]
                c2 <- mtx[g2, ]
                coexp.rates <- c(
                    coexp.rates,
                    sum(c1 > 0 & c2 > 0) / sum(c1 > 0 | c2 > 0))
            }
        }
    }
    return(coexp.rates)
}

# loop over objects and compute coexpression rates
rates <- c()
for (n in names(objects)) {
    coexp.rates <- mean(get.coexpression.rate(objects[[n]], ref.markers.ss))
    rates <- c(rates, coexp.rates)
}

# construct df containing MECR values for various gene pairs for each tech
df <- data.frame(
    Rate = rates,
    Technology = c(
        "Vizgen",
        "10x",
        "Resolve",
        "EELFISH",
        "MERFISH"),
    Fraction = c(
        sum(vizgen.counts) / vizgen.nmols,
        sum(tenx.counts) / tenx.nmols,
        sum(resolve.counts) / resolve.nmols,
        sum(eelfish.counts) / eelfish.nmols,
        sum(merfish.counts) / merfish.nmols))

# construct boxplot displaying MECR for each technology
p <- ggplot(df, aes(x = Rate, y = Fraction, color = Technology)) +
  geom_point(size = 8) + theme_bw() +
  scale_color_manual(values = colors.mapping) + NoLegend() +
  theme(text = element_text(size = 20), axis.title = element_blank()) +
  ggtitle("MECR vs. Fraction of molecules in cells")

ggsave(
    args[length(args)], p, width = 2800,
    height = 2800, units = "px", dpi = 400)
