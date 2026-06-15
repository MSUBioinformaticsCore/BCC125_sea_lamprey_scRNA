library(SingleCellExperiment)
library(ggplot2)
library(dplyr)

plot_expr_proportion <- function(sce,
                                 gene,
                                 celltype_col = "celltype",
                                 batch_col    = "batch",
                                 assay_name   = "counts",
                                 threshold    = 0,
                                 ncol         = 3) {
  
  stopifnot(gene %in% rownames(sce))
  stopifnot(celltype_col %in% names(colData(sce)))
  stopifnot(batch_col    %in% names(colData(sce)))
  
  expr <- assay(sce, assay_name)[gene, ]
  
  df <- data.frame(
    expressed = expr > threshold,
    celltype  = colData(sce)[[celltype_col]],
    batch     = colData(sce)[[batch_col]]
  )
  
  prop_df <- df %>%
    group_by(batch, celltype) %>%
    summarise(
      proportion = mean(expressed),
      n_cells    = n(),
      .groups    = "drop"
    )
  
  ggplot(prop_df, aes(x = batch, y = proportion, fill = batch)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(
      aes(label = scales::percent(proportion, accuracy = 0.1)),
      vjust = -0.4, size = 3
    ) +
    facet_wrap(~ celltype, ncol = ncol) +
    scale_y_continuous(
      labels = scales::percent_format(),
      limits = c(0, 1.20),
      breaks = seq(0, 1, by = 0.25)
    ) +
    labs(
      title = paste("Proportion of cells expressing", gene),
      x     = "Batch",
      y     = "% cells expressing",
      fill  = "Batch"
    ) +
    theme_bw() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}
plot_expr_proportion(combined.sce,
                     gene         = "LOC116951502",
                     celltype_col = "CellType",
                     batch_col    = "batch",
                     assay_name   = "logcounts",
                     threshold    = 0,
                     ncol=2)