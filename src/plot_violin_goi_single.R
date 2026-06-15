plot_violin_goi_single <- function(sce, markers, goi, gene_description, anno, palette, title = NULL) {
  
  fdr_to_stars <- function(fdr) {
    dplyr::case_when(
      fdr < 0.001 ~ "***",
      fdr < 0.01  ~ "**",
      fdr < 0.05  ~ "*",
      TRUE        ~ ""
    )
  }
  
  # extract and remap plot data
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
      X = factor(X, levels = rev(anno))  # rev so top of legend = top of plot
    )
  
  # y positions per cluster
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
    coord_flip(clip = "off", ylim = c(0, NA)) +
    labs(x = NULL, y = "Log-normalized expression", title = title) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x      = element_text(color = "#3a3a3a"),
      axis.ticks.x     = element_blank(),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.title       = element_text(color = "#2e2e2e"),
      axis.line        = element_line(color = "#aaaaaa", linewidth = 0.4),
      axis.ticks       = element_line(color = "#aaaaaa", linewidth = 0.4),
      legend.key.size  = unit(0.4, "cm")
    ) + 
    guides(fill = guide_legend(reverse = TRUE))
}

plot_violin_goi_single(
  sce            = combined.sce,
  markers        = scatter_cluster_marker_genes,
  goi            = "DDX4",
  gene_description = gene_description,
  anno           = anno,
  palette        = palette_20,
  title          = "DDX4"
)

ggsave(file = "DDX4.png", width = 5, height = 4)
