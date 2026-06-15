plot_violin_goi <- function(sce, markers, goi, gene_description, anno, palette, title = NULL) {
  
  fdr_to_stars <- function(fdr) {
    dplyr::case_when(
      fdr < 0.001 ~ "***",
      fdr < 0.01  ~ "**",
      fdr < 0.05  ~ "*",
      TRUE        ~ ""
    )
  }
  
  # Gene ID mapping
  goi_loc <- gene_description %>%
    filter(WangGeneID %in% goi) %>%
    pull(Gene)
  
  expr_mat <- as.matrix(logcounts(sce)[goi_loc, , drop = FALSE])
  rownames(expr_mat) <- goi
  
  # extract and remap plot data FIRST
  plot_df <- plotExpression(
    sce,
    features      = goi,
    scattermore   = TRUE,
    point_size    = 0,
    swap_rownames = "WangGeneID",
    x             = "label",
    ncol          = 1
  )$data %>%
    mutate(
      X = anno[as.numeric(gsub("cluster_", "", as.character(X)))],
      X = factor(X, levels = anno)
    )
  
  # y positions per Feature + cluster
  y_max_df <- plot_df %>%
    group_by(Feature, X) %>%
    summarise(y_pos = max(Y) * 1.1, .groups = "drop")
  
  expand_df <- plot_df %>%
    group_by(Feature) %>%
    summarise(y_expand = max(Y) * 1.4, .groups = "drop") %>%
    mutate(X = anno[1])
  
  # significance df
  sig_df <- markers %>%
    filter(WangGeneID %in% goi) %>%
    mutate(
      stars   = fdr_to_stars(FDR),
      Feature = WangGeneID,
      X       = anno[as.numeric(gsub("cluster_", "", as.character(cluster)))]
    ) %>%
    filter(stars != "") %>%
    left_join(y_max_df, by = c("Feature", "X"))
  
  if (length(goi) > 1) {
    gene_order <- plot_df %>%
      group_by(Feature, X) %>%
      summarise(mean_expr = mean(Y), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = X, values_from = mean_expr, values_fill = 0) %>%
      tibble::column_to_rownames("Feature") %>%
      dist() %>%
      hclust(method = "ward.D2") %>%
      { .$labels[.$order] }
    
    plot_df   <- plot_df   %>% mutate(Feature = factor(Feature, levels = gene_order))
    sig_df    <- sig_df    %>% mutate(Feature = factor(Feature, levels = gene_order))
    expand_df <- expand_df %>% mutate(Feature = factor(Feature, levels = gene_order))
  }
  
  cluster_colors <- setNames(palette[seq_along(anno)], anno)
  
  ggplot(plot_df, aes(x = X, y = Y)) +
    geom_violin(
      aes(fill = X, color = after_scale(fill)),
      linewidth = 0.3,
      scale     = "width",
      trim      = TRUE
    ) +
    geom_boxplot(
      width         = 0.12,
      fill          = "white",
      color         = "grey40",
      linewidth     = 0.3,
      outlier.shape = NA
    ) +
    geom_blank(
      data        = expand_df,
      aes(x = X, y = y_expand),
      inherit.aes = FALSE
    ) +
    geom_text(
      data        = sig_df,
      aes(x = X, y = y_pos, label = stars),
      inherit.aes = FALSE,
      size        = 2,
      vjust       = 0,
      fontface    = "bold",
      color       = "#2e2e2e"
    ) +
    scale_fill_manual(values = cluster_colors, name = "Cell Type") +
    coord_cartesian(clip = "off", ylim = c(0, NA)) +
    facet_wrap(~ Feature, ncol = 1, scales = "free_y") +
    labs(x = NULL, y = "Log-normalized expression", title = title) +
    theme_classic(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "#dce8ef", color = NA),
      strip.text       = element_text(color = "#2e4a5a", face = "bold", size = 10),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.text.y      = element_text(color = "#3a3a3a"),
      axis.title       = element_text(color = "#2e2e2e"),
      axis.line        = element_line(color = "#aaaaaa", linewidth = 0.4),
      axis.ticks       = element_line(color = "#aaaaaa", linewidth = 0.4),
      panel.spacing    = unit(0.8, "lines"),
      legend.key.size  = unit(0.4, "cm")
    )
}

plot_violin_goi(
  sce            = combined.sce,
  markers        = scatter_cluster_marker_genes,
  goi            = goi,
  gene_description = gene_description,
  anno           = anno,
  palette        = palette_20,
  title          = "Germ markers"
)



anno = c( 
  "1. Fibroblasts (migrating)",
  "2. Mesonephric cells",
  "3. Granulosa cells (migrating)",
  "4. Granulosa cells",
  "5. Oogonia/Oocytes",
  "6. Endothelial cells",
  "7. Primordial germ cells (migrating)",
  "8. Erythrocytes",
  "9. Pre-granulosa cells (migrating)",
  "10. Leukocytes",
  "11. Hemopoietic stem cells",
  "12. Podocytes",
  "13. Macrophages",
  "14. Pre-granulosa cells (mitotic)",
  "15. Immune cells (migrating, mitotic)",
  "16. Spermatocytes (cytokinetic)",
  "17. Male germ cells (migrating)",
  "18. Spermatocytes (apoptotic)",
  "19. Spermatocytes (GPCR-activated)"
)

plot_violin_goi(
  sce            = combined.sce,
  markers        = scatter_cluster_marker_genes,
  goi            = goi,
  gene_description = gene_description,
  anno           = anno,
  palette        = palette_20,
  title          = ""
)

ggsave(file = "HYKK_all_cells.png", width = 6, height = 8)

FOXL2 
SOX9


