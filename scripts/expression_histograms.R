# generate histograms of gene panel expression levels

library(Seurat)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

lin <- readRDS(args[1])
obj <- readRDS(args[2])

length(rownames(obj))
length(rownames(lin))

shared.genes <- intersect(rownames(obj), rownames(lin))

mtx <- lin@assays$RNA@counts
sums <- rowSums(mtx)

# Calculate the quantile boundaries
quantiles <- quantile(sums, probs=seq(0, 1, length.out = 51))
quantiles = quantiles + seq_along(quantiles) * .Machine$double.eps

# Use the cut function to assign each value to a quantile
quantile_labels <- cut(sums, breaks=quantiles, labels=1:50, include.lowest=TRUE)

# Combine the values and their respective quantiles into a dataframe
df <- data.frame(Value=sums, Quantile=as.numeric(as.character(quantile_labels)))
df$gene <- rownames(df)

# Plot a histogram of the quantiles
p <- ggplot(df[df$gene %in% shared.genes, ], aes(x=Quantile)) +
  geom_histogram(binwidth=1, fill="blue", color="white", alpha=0.7) +
  labs(title="Histogram of Quantiles", x="Quantile", y="Count") + theme_bw() + theme(text = element_text(size=30, color="black")) + ggtitle("")

ggsave(args[length(args)-1], p, units = "px",  width = 3000, height = 2800)

# avg expression
avg.exp <- sums / ncol(mtx)
df <- data.frame(gene = shared.genes, mean = log1p(avg.exp[shared.genes]))
print(df)
p <- ggplot(df, aes(x=mean)) + geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +
  theme_bw() + theme(text = element_text(size=30, color="black")) + ggtitle("") + xlim(0, 5)
ggsave(args[length(args)], p, units = "px",  width = 3000, height = 2800)
