library(ggplot2)
library(dplyr)
library(tidyr)

plot_pseudotime_paths <- function(sce, pseudo_mat, goi, gene_description,
                                  tradeseq_res,
                                  paths = c("Oocytes", "SPC.cytokinetic"),
                                  title = NULL) {
  
  # Resolve gene IDs
  goi_loc <- gene_description %>%
    filter(WangGeneID %in% goi) %>%
    pull(Gene)
  
  expr_mat <- as.matrix(logcounts(sce)[goi_loc, , drop = FALSE])
  rownames(expr_mat) <- gene_description %>%
    filter(Gene %in% goi_loc) %>%
    select(Gene, WangGeneID) %>%
    deframe() %>%
    { .[goi_loc] }
  
  # Pseudotime long format
  pseudo_long <- pseudo_mat[, paths, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("Cell") %>%
    tidyr::pivot_longer(-Cell, names_to = "Path", values_to = "Pseudotime") %>%
    filter(!is.na(Pseudotime))
  
  # Expression long format
  expr_long <- as.data.frame(t(expr_mat)) %>%
    rownames_to_column("Cell") %>%
    tidyr::pivot_longer(-Cell, names_to = "Gene", values_to = "Expression")
  
  # Join
  plot_df <- pseudo_long %>%
    left_join(expr_long, by = "Cell") %>%
    filter(!is.na(Expression))
  
  # FDR labels — one per Gene + Path combination
  fdr_labels <- tradeseq_res %>%
    filter(WangGeneID %in% goi, Path %in% paths) %>%
    mutate(
      label = paste0(Path, ": ", formatC(FDR, format = "e", digits = 2))
    ) %>%
    group_by(WangGeneID) %>%
    summarise(fdr_text = paste(label, collapse = "\n"), .groups = "drop") %>%
    dplyr::rename(Gene = WangGeneID)
  
  # Path colors
  path_colors <- setNames(
    c("#E377C2", "#1F77B4"),
    paths
  )
  
  ggplot(plot_df, aes(x = Pseudotime, y = Expression, color = Path, fill = Path)) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                se = TRUE, alpha = 0.15, linewidth = 0.8) +
    geom_text(
      data        = fdr_labels,
      aes(x = Inf, y = Inf, label = fdr_text),
      inherit.aes = FALSE,
      hjust       = 1.05,
      vjust       = 1.2,
      size        = 2.5,
      color       = "grey30",
      lineheight  = 1.4
    ) +
    scale_color_manual(values = path_colors) +
    scale_fill_manual(values = path_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.40))) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
    labs(x = "Pseudotime", y = "Log-normalized expression",
         color = "Path", fill = "Path", title = title) +
    theme_classic(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "#dce8ef", color = NA),
      strip.text       = element_text(color = "#2e4a5a", face = "bold", size = 9),
      axis.text        = element_text(color = "#3a3a3a"),
      axis.line        = element_line(color = "#aaaaaa", linewidth = 0.4),
      axis.ticks       = element_line(color = "#aaaaaa", linewidth = 0.4),
      legend.position  = "top",
      panel.spacing    = unit(0.8, "lines")
    ) 
}

plot_pseudotime_paths(
  sce              = combined.sex.sce,
  pseudo_mat       = pseudo_mat,
  goi              = goi,
  gene_description = gene_description,
  paths            = c("Oocytes", "SPC.cytokinetic"),
  title            = "Wnt5 cell pseudotime",
  tradeseq_res     = no_support_sex_pseudotime_de_genes
)

ggsave(file = "Wnt5_pseudotime.png", width = 5, height = 6)