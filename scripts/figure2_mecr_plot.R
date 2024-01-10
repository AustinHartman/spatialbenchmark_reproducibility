### The purpose of this script is to score the coexpression using exclusive markers
# from the linnarson reference
# Identify markers which are > 25% expressed in marker cell type and < 1% expression not in that cell type
# - Intersect this list with the panel of markers for Vizgen, 10x, and Resolve.
# - For each pair of non-coexpressed exclusive markers, compute the fraction of coexpressing cells.
# - Plot fraction of coexpressing cells as a histogram.
# - I would expect 10x to be the highest, then Vizgen, then Resolve.

library(Seurat)
library(ggplot2)
library(SeuratWrappers)
source("./R/utils.R")

args <- commandArgs(trailingOnly = TRUE)

# Load in other objects for comparison (not-necessarily spatial)
spatial.objects <- list()
l <- (length(args) - 2) / 2
ps <- 2
ns <- ps + l + 1
print(args)
for (a in 0:(l-1)) {
    print(a)
    print(args[ps + a])
    print(file.exists(args[ps + a]))
    print(args[ns + a])
    spatial.objects[[args[ns + a]]] <- readRDS(args[ps + a])
}

# Load reference and identify markers
ref <- readRDS(args[1])
ref <- NormalizeData(ref)
Idents(ref) <- "Class"
ref.markers <- RunPrestoAll(ref)
ref.markers.ss <- ref.markers[
  ref.markers$pct.1 > 0.25 & ref.markers$pct.2 < 0.01,
]

get.coexpression.rate <- function(obj, sc.markers) {
  coexp.rates <- c()
  genes <- intersect(rownames(obj), rownames(sc.markers))
  print(length(genes))
  if (length(genes) > 25) { genes <- sample(genes, 25) }
  mtx <- as.matrix(obj@assays$RNA@counts[genes, ])
  for (g1 in genes) {
    for (g2 in genes) {
      if ((g1 != g2) && (g1 > g2) && (sc.markers[g1, "cluster"] != sc.markers[g2, "cluster"])) {
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

rates <- c()
technology <- c()
for (n in names(spatial.objects)) {
    coexp.rates <- get.coexpression.rate(spatial.objects[[n]], ref.markers.ss)
    rates <- c(rates, coexp.rates)
    technology <- c(technology, rep(n, length(coexp.rates)))
    message(paste("Coexpression rate", n, "done"))
}
coexp.rates <- get.coexpression.rate(ref, ref.markers.ss)
rates <- c(rates, coexp.rates)
technology <- c(technology, rep("Linnarsson", length(coexp.rates)))
df <- data.frame(Rates = rates, Technology = technology)
technologies <- unique(df$Technology)
if ("Baysor" %in% technologies) {
  colors.mapping <- list(
    "Nuclear" = "#FF5733",
    "Cell" = "#59C230",
    "Baysor" = "#3258A5",
    "Linnarsson" = "black",
    "Allen" = "black")
} else {
  colors.mapping <- colors[tolower(technologies)]
  colors.mapping <- lapply(colors.mapping, function(x) { ifelse(is.null(x), "black", x) })
  names(colors.mapping) <- technologies
}
plot <- ggplot(df, aes(x = Rates, y = Technology, fill = Technology)) +
  geom_boxplot() + theme_bw() +
  scale_fill_manual(values = colors.mapping) + NoLegend() +
  theme(text = element_text(size = 25), axis.title = element_blank()) +
  ggtitle("Coexpression rate")
message(paste("Saving plot to", args[length(args)]))
ggsave(args[length(args)], plot, width = 2800, height = 2800, units = "px", dpi = 400)
