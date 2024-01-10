library(Seurat)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
message(paste0("Arguments:", args))

object <- readRDS(file.path(getwd(), args[1]))
object2 <- readRDS(file.path(getwd(), args[2]))

df <- data.frame(
  rep1 = log10(rowSums(object@assays$RNA@counts)),
  rep2 = log10(rowSums(object2@assays$RNA@counts))
)

min <- min(min(df$rep1), min(df$rep2))
max <- max(max(df$rep1), max(df$rep2))
t <- paste0(args[3], " replicate gene counts")
p <- ggplot(df, mapping = aes(rep1, rep2)) +
  geom_point() + theme_classic() +
  xlab("Replicate 1 - log10 scale") + ylab("Replicate 2 - log10 scale") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme(
    text = element_text(size=25, color="black")
  ) + xlim(min, max) + ylim(min, max) + ggtitle(t)

print(paste("Correlation", args[3], cor(df$rep1, df$rep2)))

save.path <- file.path("./data/", paste0("replicate_plot_", args[3], ".png"))
ggsave(filename = save.path, plot = p, units = "px", width = 2500, height = 2500)
