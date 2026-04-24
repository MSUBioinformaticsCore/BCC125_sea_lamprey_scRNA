markers_wide = read.csv(paste0(data_dir, "/Cell_Type_Marker_Genes_12-9-25.csv"))
markers_wide = markers_wide[,2:ncol(markers_wide)]
  
markers_long <- markers_wide %>%
  pivot_longer(
    cols = everything(),
    names_to = "CellType",
    values_to = "Gene",
    values_drop_na = TRUE
  ) %>%  mutate(
    main_part = str_remove(Gene, "\\s*\\(.*\\)"),  # Remove parentheses and content
    paren_part = str_extract(Gene, "(?<=\\().*(?=\\))")  # Extract content inside parentheses
  ) %>%
  separate_rows(main_part, sep = "/") %>%
  mutate(main_part = trimws(main_part),
         paren_part = trimws(paren_part))

marker_df = markers_long %>%
  select(CellType, main_part) %>%
  dplyr::rename(WangGeneID = "main_part") %>%
  filter(!WangGeneID == "") %>%
  left_join(gene_description) %>%
  filter(!is.na(Gene))

library(tibble)
# -------------------------------------------
# 1. Expression matrix scaled
# -------------------------------------------

# germ cells
sce = combined.sex.sce
colData(sce) = colData(merged.sex.sce)
markers = read.csv(paste0(results_dir, "/no.support_sex_scatter_cluster_marker_genes.csv"), row.names = 1)

anno =
  c("1. Epithelial cells I",
    "2. Epithelial cells II",
    "3. Mitotic spermatogonia",
    "4. Oocytes",
    "5. Epithelial cells III",
    "6. Mitotic oogonia",
    "7. Mitotic PGCs",
    "8. Spermatogonia",
    "9. Mesonephric cells I",
    "10. Spermatocytes I",
    "11. Mesonephric cells II",
    "12. Male germ cells",
    "13. Spermatocytes I",
    "14. Spermatocytes II",
    "15. Migrating PGCs"
  )

# # all cells
sce = combined.sce
markers = read.csv("/mnt/ufs18/rs-013/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA/results/20250930/scatter_cluster_marker_genes.csv", row.names = 1)

anno = c( "1. Fibroblasts",
          "2. Mesonephric cells",
          "3. Granulosa cells I",
          "4. Granulosa cells II",
          "5. PGCs/Oogonia/Oocytes",
          "6. Endothelial cells",
          "7. Migrating PGCs",
          "8. Erythrocytes",
          "9. Pre-granulosa cells I",
          "10. Leukocytes",
          "11. Hemopoietic stem cells",
          "12. Podocytes",
          "13. Macrophages",
          "14. Pre-granulosa cells II",
          "15. Migrating immune cells",
          "16. Spermatocytes I",
          "17. Spermatocytes II",
          "18. Spermatogonia I",
          "19. Spermatogonia II")

colData(sce) = colData(merged.sce)

ct = anno[colData(sce)$label]
sce$CellType = ct
sce$CellType = factor(sce$CellType, levels= anno)
sce$cluster = sce$label

cluster_index = as.numeric(gsub("cluster_","",markers$cluster))
markers$CellType = anno[cluster_index]

markers = 
  markers %>%
  select(WangGeneID, CellType, FDR, summary.logFC) %>%
  dplyr::rename(Gene = "WangGeneID")

known_markers = chung_davidson_markers

nested_df <- known_markers %>%
  group_by(CellType) %>%
  tidyr::nest()

# Then extract as a named list
list_of_dfs <- nested_df$data
names(list_of_dfs) <- nested_df$CellType
names(list_of_dfs) = gsub("\\(", "", names(list_of_dfs))
file_names = gsub(" ", "_", names(list_of_dfs))

for(i in 1:length(list_of_dfs)){

print(names(list_of_dfs)[i])
goi = list_of_dfs[[i]]$Gene
goi_in = goi[goi %in% rownames(sce)]
goi_notin = goi[!goi %in% rownames(sce)]

notin_id = 
  gene_description %>%
  filter(WangGeneID %in% goi_notin) %>%
  pull(Gene)

goi = c(goi_in, notin_id)

sce.plot = sce[goi,]

goi_names = 
  gene_description %>%
  filter(Gene %in% goi) %>%
  select(Gene, WangGeneID) %>%
  deframe()

goi_names = goi_names[goi]

rownames(sce.plot) = goi_names

mat <- logcounts(sce.plot)
mat = as.matrix(mat)

meta <- as.data.frame(colData(sce.plot))
meta <- meta %>% arrange(cluster)

df_long <- as.data.frame(mat) %>%
  rownames_to_column("Gene") %>%
  tidyr::pivot_longer(-Gene, names_to = "Cell", values_to = "Expression") %>%
  left_join(meta %>% tibble::rownames_to_column("Cell"), by = "Cell") %>%
  select(Gene, Cell, CellType, Expression)

plot_data <- df_long %>%
  group_by(Gene, CellType) %>%
  summarize(
    percent_expressed = mean(Expression > 0) * 100,
    avg_expression = mean(Expression[Expression > 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Gene) %>%
  mutate(
    scaled_expression = scale(avg_expression)[,1]  # Z-score by gene
  ) %>%
  ungroup() %>%
  mutate(
    text_color = ifelse(avg_expression > 5,
                        "black", "white")
  ) %>%
  left_join(markers)

percent_cells = 
  df_long %>%
  group_by(Gene, CellType) %>%
  summarize(
    percent_expressing = mean(Expression > 0) * 100,
    n_cells = n(),
    n_expressing = sum(Expression > 0),
    .groups = "drop"
  )

percent_cells_wide = 
  percent_cells %>%
  tidyr::pivot_wider(id_cols = "Gene", 
              names_from="CellType", 
              values_from = "percent_expressing")

gene_rows = percent_cells_wide$Gene
percent_cells_wide$Gene = NULL
percent_cells_mat = as.matrix(percent_cells_wide)
rownames(percent_cells_mat) = gene_rows 

if(length(goi) > 2){
  gene_dist <- dist(percent_cells_mat)
  gene_clust <- hclust(gene_dist, method = "ward.D2")
  gene_order <- gene_clust$labels[gene_clust$order]
  percent_cells$Gene = factor(percent_cells$Gene, levels = gene_order)
  plot_data$Gene = factor(plot_data$Gene, levels = gene_order)
  
}

# -------------------------------------------
# Calculate dynamic height
# -------------------------------------------
n_genes <- length(unique(plot_data$Gene))
max_bubble_size <- 5  # from scale_size_continuous range
bubble_size_mm <- max_bubble_size * 0.35278  # convert points to mm (approximate)
gene_height <- bubble_size_mm / 25.4  # convert mm to inches
plot_height <- max(n_genes * gene_height * 2 + 1, 7)  # minimum height of 8 inches

# -------------------------------------------
# 6. Colors
# -------------------------------------------

all_label_colors = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78",
                     "#2CA02C", "#98DF8A", "#D62728", "#FF9896", 
                     "#9467BD", "#C5B0D5", "#8C564B", "#C49C94", 
                     "#E377C2", "#F7B6D2", "#7F7F7F", "#C7C7C7", 
                     "#BCBD22", "#DBDB8D", "#17BECF")

label_colors = all_label_colors[1:length(unique(percent_cells$CellType))] 

names(label_colors) = unique(percent_cells$CellType)

# -------------------------------------------
# 7. Annotation bar
# -------------------------------------------

anno_bar <- ggplot(percent_cells, aes(x = CellType, y = 1, fill = CellType)) +
  geom_tile() +
  scale_fill_manual(values = c(label_colors, "NA" = "white"), na.value = "white", drop = FALSE) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_minimal(base_size = 10) +
  theme(
    #legend.position = "none",
    legend.key.size = unit(.4, "cm"),
    legend.key.width = unit(.4, "cm"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(0, 5, 0, 40)
  ) + labs(title = names(list_of_dfs)[i])

# -------------------------------------------
# 8. Heatmap
# -------------------------------------------
plot_data$CellType = factor(plot_data$CellType, levels = anno)

library(ggnewscale)

bubble = 
  ggplot(plot_data, aes(x = CellType, y = Gene)) +
  geom_point(aes(size = percent_expressed, color = avg_expression)) +
  scale_size_continuous(
    name = "% Expressed",
    range = c(1, 5)  # min and max bubble size
  ) +
  scale_color_viridis_c(option = "C",
                        name = "Avg Expression") +
  new_scale_color() + 
  geom_text(
    data = plot_data %>% filter(FDR<.01, summary.logFC > 0.5),  
    aes(label = "*", color = text_color),
    size = 2,
    vjust = 0.6,
    show.legend = FALSE
  ) +
  scale_color_identity() +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    axis.text.y = element_text(size = 7, 
                               color = "black",
                               face = "italic",
                               family = "sans"),
  )

# -------------------------------------------
# 9. Combine
# -------------------------------------------
anno_bar / bubble + plot_layout(heights = c(0.05, 1), guides = "collect")
ggsave(file = paste0(results_dir, "/20260106_known_marker_plots/", file_names[i], "_germ_cells_bubble.png"), 
       width = 8, height=plot_height, limitsize=FALSE)
}
