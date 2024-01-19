# The purpose of this plot is to generate counts per feature rank plots for
# each of the available technologies with background probes highlighted.

library(Seurat)
library(ggplot2)

color.mapping <- list(
    "Gene probes" = "black",
    "Negative probes" = "firebrick1",
    "Negative control probes" = "firebrick1",
    "Negative control codewords" = "firebrick1",
    "Blanks" = "firebrick1")

args <- commandArgs(trailingOnly = TRUE)
ratios <- list()
for (i in seq(1, 9, 2)) {
    print(args[i])
    print(args[i + 1])
    print("")
    molecules <- read.csv(args[i])
    technology <- args[i + 1]
    counts <- table(molecules$gene)

    if (technology == "Vizgen") {
        background.features <- names(counts)[grepl("Blank*", names(counts))]
    } else if (technology == "10x") {
        blank.genes <- names(counts)[grepl("BLANK*", names(counts))]
        neg.cw <- names(counts)[grepl("NegControlCodeword*", names(counts))]
        neg.con.probes <- names(counts)[
            grepl("NegControlProbe*", names(counts))]
    } else if (technology == "Resolve") {
        background.features <- names(counts)[grepl("^FP.*", names(counts))]
    } else if (technology == "EEL FISH") {
        background.features <- names(counts)[grepl("Control*", names(counts))]
    } else if (technology == "MERFISH") {
        background.features <- names(counts)[grepl("blank*", names(counts))]
    }

    df <- data.frame(
        counts = unlist(counts),
        category = "Gene probes",
        order = 1)
    colnames(df) <- c("features", "counts", "category", "order")
    if (technology == "10x") {
        df[df$features %in% blank.genes, "category"] <- "Blanks"
        df[df$features %in% blank.genes, "order"] <- 4

        df[df$features %in% neg.cw, "category"] <- "Negative control codewords"
        df[df$features %in% neg.cw, "order"] <- 3

        df[df$features %in% neg.con.probes, "category"] <- "Negative control probes"
        df[df$features %in% neg.con.probes, "order"] <- 2

        ratio <- mean(df[df$category == "Negative control probes", "counts"]) /
            mean(df[df$category == "Gene probes", "counts"])
    } else {
        df[
            df$features %in% background.features,
            "category"] <- "Negative probes"
        ratio <- mean(df[df$category == "Negative probes", "counts"]) /
            mean(df[df$category == "Gene probes", "counts"])
    }
    ratios[[technology]] <- ratio
}

df <- data.frame(
    technology = names(ratios),
    ratio = unlist(ratios))
p <- ggplot(df, aes(y = reorder(technology, ratio), x = ratio)) +
    geom_bar(stat = "identity", fill = "black") +
    coord_flip() +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        text = element_text(size = 30)) +
    labs(
        title = "Ratio of mean counts per gene to mean counts per background probe",
        subtitle = "Background probes are defined as probes that do not target any gene",
        x = "Ratio of mean counts per background probe to mean counts per gene",
        y = "Technology") +
    scale_x_continuous(expand = c(0, 0))

ggsave(args[length(args)], p, units = "px",  width = 3000, height = 3800)
