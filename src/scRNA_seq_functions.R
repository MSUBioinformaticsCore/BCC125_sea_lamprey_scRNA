# loadRData ---------------------------------------------------------------
#' @name loadRData
#' @description loads an RData file, and returns it
#' @param fileName path to RData file
# https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# phyperOverlapSets ------------------------------------------------------------
#' @name phyperOverlapSets
#' @description performs pairwise hypergeometric tests between multiple sets 
#' of genes i.e. cluster marker genes from two different single cell
#' different single cell data sets
#' @param table1 a dataframe with two columns, 
#' ie gene and cluster assignment, from dataset 1      
#' @param table2 a dataframe with two columns, 
#' ie gene and cluster assignment, from dataset 2
#' format example
#' --------------  
#' Gene  Cluster
#' Gene1 cluster1
#' Gene2 cluster1
#' Gene3 cluster2
#' Gene4 cluster3
#' @param background the number of genes in the gene universe 
#' ie. all expressed genes
#' @return a list with two elements: 
#' [1]intersecting.genes: a tidy table where each row shows two groups 
#' (one from table1 and one from table2) and an intersecting gene or "none"
#' [2] pval.scores: a tidy table where each row is a pair of groups 
#' (one from table1 and one from table2) with the number of genes in each group,
#' the nunber of intersecting genes, the uncorrected pval, and the enrichment score

phyperOverlapSets <- function(table1, table2, background){
  require(dplyr)
  
  # make table1 into a list divided by group
  colnames(table1) = c("Gene","Group")
  tab1_groups = unique(table1$Group)
  tab1_list = list()
  
  for(group in tab1_groups){
    tab1_list[[group]] = table1 %>%
      filter(Group == group)
  }
  
  # make table2 into a list divided by group
  colnames(table2) = c("Gene","Group")
  tab2_groups = unique(table2$Group)
  tab2_list = list()
  
  for(group in tab2_groups){
    tab2_list[[group]] = table2 %>%
      filter(Group == group)
  }
  
  # make a list of dfs with the names of genes shared between
  # every group in table1 and every group in table2
  INT = list()
  # make a list of dfs with the number of genes shared between
  # every group in table1 and every group in table2
  nINT = list()
  
  for(group in tab1_groups){
    genes = tab1_list[[group]]$Gene
    INT_list = lapply(tab2_list, 
                      function(x){
                        y = intersect(x$Gene, genes);
                        if(length(y) > 0){
                          z = data.frame(SharedGene = y)}
                        else{z = data.frame(SharedGene = "none")};
                        z
                      })
    
    for(i in 1:length(INT_list)){
      INT_list[[i]]$Group2 = names(INT_list)[i]
    }
    
    INT[[group]] = do.call(rbind,INT_list)
    INT[[group]]$Group1 = group
    
    INT[[group]] = INT[[group]] %>%
      select(Group1, Group2, SharedGene)
    
    nINT_list = lapply(tab2_list, 
                       function(x){
                         y = intersect(x$Gene, genes);
                         z = length(y);
                         z
                       })
    
    nINT_vec = unlist(nINT_list)
    nGroup2 = unlist(lapply(tab2_list, nrow))
    
    nINT[[group]] = data.frame(Group1 = group,
                               nGroup1 = length(genes),
                               Group2 = names(nINT_vec),
                               nGroup2  = nGroup2,
                               nSharedGenes = nINT_vec)
    
    
  }
  
  INTdf = do.call(rbind, INT)
  nINTdf = do.call(rbind, nINT)
  
  nINTdf = nINTdf %>%
    mutate(Pval = phyper(nSharedGenes-1,
                         nGroup2,
                         background-nGroup2,
                         nGroup1,
                         lower.tail= FALSE)
    )
  
  nINTdf = nINTdf %>%
    mutate(Enrichment = 
             log2(
               (nSharedGenes/nGroup2)/
                 (nGroup1/background))
    )
  
  return(list(intersecting.genes = INTdf, 
              pval.scores = nINTdf))
}

# fgseaEnrichment ---------------------------------------------------------
#' @name fgseaEnrichment
#' @description performs GSEA enrichments pairwise between a ranked lists of  
#' of genes i.e. cluster marker genes with some stat like -log10FDR and other gene sets
#' @param stat_table a dataframe with three columns from dataset 1, 
#' ie Gene, Set, Stat. Should include all genes with no filter.
#' format example
#' --------------  
#' Gene  Set      Stat
#' Gene1 cluster1 10
#' Gene2 cluster1 100
#' Gene3 cluster2 5
#' Gene4 cluster3 20
#' @param set_table a dataframe with two columns, 
#' ie gene and cluster assignment, from dataset 2
#' format example
#' --------------  
#' Gene  Set
#' Gene1 cluster1
#' Gene2 cluster1
#' Gene3 cluster2
#' Gene4 cluster3
#' @param scoreType "pos" if all stats are positive not signed, "std" otherwise
#' ie. all expressed genes
#' @return  a tidy table where each row is a pair of groups 
#' (one from stat_table and one from set_table) with the number of genes in each group,
#' the number of intersecting genes, the uncorrected pval, and the enrichment score

fgseaEnrichment = function(stat_table, set_table, scoreType) {
  
  require(fgsea)
  require(tidyverse)
  
  set_list = list()
  for(set in unique(set_table$Set)){
    
    set_list[[set]] = 
      set_table %>%
      filter(Set == set) %>%
      pull(Gene)
  } 
  
  stat_list = list()
  for(set in unique(as.character(stat_table$Set))){
    
    stat_list[[set]] = 
      stat_table %>%
      filter(Set == set) %>%
      pull(Stat)
    
    stat_list[[set]] =  
      as.numeric(stat_list[[set]])
    
    max_stat = max(stat_list[[set]][is.finite(stat_list[[set]])])
    stat_list[[set]][is.infinite(stat_list[[set]])] = max_stat + 1
    
    names(stat_list[[set]]) =
      stat_table %>%
      filter(Set == set) %>%
      pull(Gene)
  }
  
  apply_fgseaMultilevel = function(stats, pathways, minSize, scoreType){
    
    res = 
      fgseaMultilevel(
        pathways = pathways,
        stats = stats,
        minSize = minSize,
        scoreType = scoreType
      )
    
    return(res)
  }
  
  res_list = lapply(stat_list, 
                    apply_fgseaMultilevel,
                    pathways = set_list,
                    minSize = 1,
                    scoreType = scoreType)
  
  names(res_list) = names(stat_list)
  
  for(i in 1:length(res_list)){
    
    res_list[[i]]$StatDataset = names(res_list)[i]
    
    res_list[[i]] = 
      res_list[[i]] %>%
      select(-leadingEdge) %>%
      rename(SetDataset = pathway) %>%
      select(StatDataset, everything())
  }
  
  res = do.call(rbind, res_list)
  return(res)
}

# find_stable_k -----------------------------------------------------------

find_stable_k <- function(sce, k_range = seq(5, 50, by = 5), 
                          n_iterations = 10, use.dimred = "PCA", seed = 42) {
  
  set.seed(seed)
  
  # Store results
  results <- list()
  
  for (k in k_range) {
    cat("Testing k =", k, "\n")
    
    # Store cluster assignments across iterations
    cluster_matrix <- matrix(0, nrow = ncol(sce), ncol = n_iterations)
    n_clusters <- numeric(n_iterations)
    
    for (i in 1:n_iterations) {
      # Build SNN graph
      snn <- buildSNNGraph(sce, k = k, use.dimred = use.dimred)
      
      # Louvain clustering
      clusters <- cluster_louvain(snn)$membership
      cluster_matrix[, i] <- clusters
      n_clusters[i] <- length(unique(clusters))
    }
    
    # Calculate stability metrics
    # 1. Variation in number of clusters
    cluster_count_sd <- sd(n_clusters)
    cluster_count_mean <- mean(n_clusters)
    cluster_count_cv <- cluster_count_sd / cluster_count_mean
    
    # 2. Average pairwise ARI (Adjusted Rand Index) between iterations
    ari_values <- numeric()
    for (i in 1:(n_iterations - 1)) {
      for (j in (i + 1):n_iterations) {
        ari <- mclust::adjustedRandIndex(cluster_matrix[, i], 
                                         cluster_matrix[, j])
        ari_values <- c(ari_values, ari)
      }
    }
    mean_ari <- mean(ari_values)
    
    results[[as.character(k)]] <- list(
      k = k,
      mean_n_clusters = cluster_count_mean,
      sd_n_clusters = cluster_count_sd,
      cv_n_clusters = cluster_count_cv,
      mean_ari = mean_ari,
      min_ari = min(ari_values),
      cluster_matrix = cluster_matrix
    )
  }
  
  return(results)
}

# plot_stability_metrics --------------------------------------------------
# Function to visualize stability results
plot_stability_metrics <- function(results) {
  
  # Extract metrics into data frame
  df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      k = x$k,
      mean_n_clusters = x$mean_n_clusters,
      cv_n_clusters = x$cv_n_clusters,
      mean_ari = x$mean_ari,
      min_ari = x$min_ari
    )
  }))
  
  # Create plots
  p1 <- ggplot(df, aes(x = k, y = mean_n_clusters)) +
    geom_line(color = "steelblue", size = 1) +
    geom_point(size = 3) +
    labs(title = "Number of Clusters vs k",
         x = "k (number of neighbors)",
         y = "Mean number of clusters") +
    theme_minimal()
  
  p2 <- ggplot(df, aes(x = k, y = mean_ari)) +
    geom_line(color = "darkgreen", size = 1) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
    labs(title = "Clustering Stability (ARI) vs k",
         x = "k (number of neighbors)",
         y = "Mean ARI between iterations",
         subtitle = "Higher ARI = more stable (aim for > 0.9)") +
    ylim(0, 1) +
    theme_minimal()
  
  p3 <- ggplot(df, aes(x = k, y = cv_n_clusters)) +
    geom_line(color = "darkorange", size = 1) +
    geom_point(size = 3) +
    labs(title = "Coefficient of Variation vs k",
         x = "k (number of neighbors)",
         y = "CV of cluster counts",
         subtitle = "Lower CV = more stable") +
    theme_minimal()
  
  return(list(p1 = p1, p2 = p2, p3 = p3, data = df))
}


# get_avg_expression ------------------------------------------------------
get_avg_expression <- function(sce, genes, cluster_col, assay = "logcounts") {
  # Filter to genes present in object
  genes <- genes[genes %in% rownames(sce)]
  if (length(genes) == 0) stop("No genes found in SCE object")
  
  expr_mat <- assay(sce, assay)[genes, , drop = FALSE]
  clusters <- colData(sce)[[cluster_col]]
  
  # Average per cluster
  avg <- sapply(unique(clusters), function(cl) {
    rowMeans(expr_mat[, clusters == cl, drop = FALSE])
  })
  return(avg)
}


# plot_pseudotime_heatmap -------------------------------------------------

plot_pseudotime_heatmap <- function(
    sce_object,
    gene_ids,                # e.g. top50
    gene_labels,             # e.g. top50_wang
    cell_subset,             # e.g. on.path8
    pseudotime_column = "Pseudotime",
    heatmap_limits = c(-4, 5),
    title = NULL 
) {
  
  # -------------------------------------------
  # 1. Expression matrix centered
  # -------------------------------------------
  mat <- logcounts(sce_object[gene_ids, cell_subset])
  rownames(mat) <- gene_labels
  
  mat_scaled <- t(scale(t(mat), center = TRUE))
  
  # -------------------------------------------
  # 2. Metadata + ordering by pseudotime
  # -------------------------------------------
  meta <- as.data.frame(colData(sce_object[, cell_subset]))
  
  # Ensure numeric ordering
  meta[[pseudotime_column]] <- as.numeric(meta[[pseudotime_column]])
  
  meta <- meta %>%
    dplyr::arrange(.data[[pseudotime_column]])
  
  # Reorder expression matrix columns
  mat_scaled <- mat_scaled[, rownames(meta)]
  
  # -------------------------------------------
  # 3. Cluster genes (rows)
  # -------------------------------------------
  gene_dist  <- dist(mat_scaled)
  gene_clust <- hclust(gene_dist, method = "ward.D2")
  gene_order <- gene_clust$labels[gene_clust$order]
  
  mat_scaled <- mat_scaled[gene_order, , drop = FALSE]
  
  # -------------------------------------------
  # 4. Long format data
  # -------------------------------------------
  df_long <- as.data.frame(mat_scaled) %>%
    tibble::rownames_to_column("Gene") %>%
    tidyr::pivot_longer(-Gene, names_to = "Cell", values_to = "Expression") %>%
    dplyr::left_join(
      meta %>% tibble::rownames_to_column("Cell"),
      by = "Cell"
    )
  
  df_long$Gene <- factor(df_long$Gene, levels = gene_order)
  
  # -------------------------------------------
  # 5. Annotation bar (pseudotime)
  # -------------------------------------------
  meta$Cell <- factor(rownames(meta), levels = rownames(meta))
  
  m_anno_bar <- ggplot2::ggplot(
    meta,
    ggplot2::aes(
      x = Cell,
      y = 1,
      fill = .data[[pseudotime_column]]
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(
      option = "viridis",
      na.value = "white",
      limits = c(0, 1)
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      legend.position = "right",
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(0, 5, 0, 40)
    ) +
    ggplot2::labs(fill = "Pseudotime") +
    labs(title = title)
  
  # -------------------------------------------
  # 6. Heatmap
  # -------------------------------------------
  df_long$Cell <- factor(df_long$Cell, levels = levels(meta$Cell))
  
  m_heatmap_plot <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = Cell, y = Gene, fill = Expression)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "steelblue3",
      mid = "#FFFFE0",
      high = "firebrick2",
      midpoint = 0,
      na.value = "white",
      limits = heatmap_limits
    ) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        size = 7,
        color = "black",
        face = "italic",
        family = "sans"
      ),
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 5, 0, 40)
    ) +
    ggplot2::labs(fill = "Scaled\nexpression")
  
  # -------------------------------------------
  # 7. Combine plots
  # -------------------------------------------
  combined_plot <- m_anno_bar / m_heatmap_plot +
    patchwork::plot_layout(
      heights = c(0.025, 1),
      guides = "collect"
    )
  
  return(combined_plot)
}


# plot_marker_bubbles -----------------------------------------------------

plot_marker_bubbles <- function(
    sce,
    markers,
    gene_description,
    cell_type_name,
    cell_type_genes,
    anno,
    results_dir,
    file_name = NULL
) {
  
  # -------------------------------------------
  # 1. Get genes of interest
  # -------------------------------------------
  goi <- cell_type_genes$Gene
  goi_in <- goi[goi %in% rownames(sce)]
  goi_notin <- goi[!goi %in% rownames(sce)]
  
  notin_id <- gene_description %>%
    filter(WangGeneID %in% goi_notin) %>%
    pull(Gene)
  
  goi <- c(goi_in, notin_id)
  
  sce.plot <- sce[goi, ]
  
  goi_names <- gene_description %>%
    filter(Gene %in% goi) %>%
    select(Gene, WangGeneID) %>%
    deframe()
  
  goi_names <- goi_names[goi]
  rownames(sce.plot) <- goi_names
  
  # -------------------------------------------
  # 2. Expression matrix + metadata
  # -------------------------------------------
  mat <- as.matrix(logcounts(sce.plot))
  meta <- as.data.frame(colData(sce.plot)) %>% arrange(cluster)
  
  df_long <- as.data.frame(mat) %>%
    rownames_to_column("Gene") %>%
    tidyr::pivot_longer(-Gene, names_to = "Cell", values_to = "Expression") %>%
    left_join(meta %>% tibble::rownames_to_column("Cell"), by = "Cell") %>%
    select(Gene, Cell, CellType, Expression)
  
  # -------------------------------------------
  # 3. Plot data
  # -------------------------------------------
  plot_data <- df_long %>%
    group_by(Gene, CellType) %>%
    summarize(
      percent_expressed = mean(Expression > 0) * 100,
      avg_expression = mean(Expression[Expression > 0], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(Gene) %>%
    mutate(scaled_expression = scale(avg_expression)[, 1]) %>%
    ungroup() %>%
    mutate(text_color = ifelse(avg_expression > 5, "black", "white")) %>%
    left_join(markers, by = c("Gene", "CellType"))
  
  # -------------------------------------------
  # 4. Percent cells + gene clustering
  # -------------------------------------------
  percent_cells <- df_long %>%
    group_by(Gene, CellType) %>%
    summarize(
      percent_expressing = mean(Expression > 0) * 100,
      n_cells = n(),
      n_expressing = sum(Expression > 0),
      .groups = "drop"
    )
  
  if (length(goi) > 2) {
    percent_cells_wide <- percent_cells %>%
      tidyr::pivot_wider(id_cols = "Gene",
                         names_from = "CellType",
                         values_from = "percent_expressing")
    
    gene_rows <- percent_cells_wide$Gene
    percent_cells_wide$Gene <- NULL
    percent_cells_mat <- as.matrix(percent_cells_wide)
    rownames(percent_cells_mat) <- gene_rows
    
    gene_dist <- dist(percent_cells_mat)
    gene_clust <- hclust(gene_dist, method = "ward.D2")
    gene_order <- gene_clust$labels[gene_clust$order]
    
    percent_cells$Gene <- factor(percent_cells$Gene, levels = gene_order)
    plot_data$Gene <- factor(plot_data$Gene, levels = gene_order)
  }
  
  # -------------------------------------------
  # 5. Dynamic height
  # -------------------------------------------
  n_genes <- length(unique(plot_data$Gene))
  bubble_size_mm <- 5 * 0.35278
  gene_height <- bubble_size_mm / 25.4
  plot_height <- max(n_genes * gene_height * 2 + 1, 7)
  
  # -------------------------------------------
  # 6. Colors
  # -------------------------------------------
  all_label_colors <- c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78",
                        "#2CA02C", "#98DF8A", "#D62728", "#FF9896",
                        "#9467BD", "#C5B0D5", "#8C564B", "#C49C94",
                        "#E377C2", "#F7B6D2", "#7F7F7F", "#C7C7C7",
                        "#BCBD22", "#DBDB8D", "#17BECF")
  
  label_colors <- all_label_colors[1:length(unique(percent_cells$CellType))]
  names(label_colors) <- unique(percent_cells$CellType)
  
  # -------------------------------------------
  # 7. Annotation bar
  # -------------------------------------------
  anno_bar <- ggplot(percent_cells, aes(x = CellType, y = 1, fill = CellType)) +
    geom_tile() +
    scale_fill_manual(values = c(label_colors, "NA" = "white"), na.value = "white", drop = FALSE) +
    scale_x_discrete(expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    theme(
      legend.key.size = unit(.4, "cm"),
      legend.key.width = unit(.4, "cm"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(0, 5, 0, 40)
    ) +
    labs(title = cell_type_name)
  
  # -------------------------------------------
  # 8. Bubble plot
  # -------------------------------------------
  plot_data$CellType <- factor(plot_data$CellType, levels = anno)
  
  bubble <- ggplot(plot_data, aes(x = CellType, y = Gene)) +
    geom_point(aes(size = percent_expressed, color = avg_expression)) +
    scale_size_continuous(name = "% Expressed", range = c(1, 5)) +
    scale_color_viridis_c(option = "C", name = "Avg Expression") +
    ggnewscale::new_scale_color() +
    geom_text(
      data = plot_data %>% filter(FDR < .01, summary.logFC > 0.5),
      aes(label = "*", color = text_color),
      size = 2, vjust = 0.6, show.legend = FALSE
    ) +
    scale_color_identity() +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      legend.key.size = unit(0.5, "cm"),
      axis.text.y = element_text(size = 7, color = "black", face = "italic", family = "sans")
    )
  
  # -------------------------------------------
  # 9. Combine + save
  # -------------------------------------------
  combined <- anno_bar / bubble +
    patchwork::plot_layout(heights = c(0.05, 1), guides = "collect")
  
  if (!is.null(results_dir)) {
    fn <- if (!is.null(file_name)) file_name else gsub(" ", "_", cell_type_name)
    ggsave(
      file = file.path(results_dir, paste0(fn, "_bubble.png")),
      plot = combined,
      width = 8,
      height = plot_height,
      limitsize = FALSE
    )
  }
  
  return(combined)
}


# run_go_enrichment -------------------------------------------------------
#' @name run_go_enrichment
#' @description Runs GO hypergeometric overlap, writes CSVs, returns scores df.
#' @param gaf filtered gene annotation file (GOID, Gene columns)
#' @param gene_description gene descriptions df with Gene, WangGeneID, description
#' @param go_names GO term name table (GOID, GOName)
#' @param test_genes two-column df: Gene and group (cluster / path)
#' @param background_n number of genes in the universe
#' @param group2_label column name for the group2 variable in output (default "Clusters")
#' @param ngroup2_label column name for nGroup2 in output (default "nClusterEnrichedGenes")
#' @param overlaps_file path for intersecting genes CSV; NULL to skip
#' @param pvals_file path for GO scores CSV; NULL to skip
#' @param pval_filter if not NULL, keep only rows with Pval < pval_filter
#' @param require_shared if TRUE, keep only rows with nSharedGenes > 0

run_go_enrichment <- function(gaf, gene_description, go_names,
                               test_genes, background_n,
                               group2_label   = "Clusters",
                               ngroup2_label  = "nClusterEnrichedGenes",
                               overlaps_file  = NULL,
                               pvals_file     = NULL,
                               pval_filter    = NULL,
                               require_shared = FALSE) {
  require(dplyr)

  go_ol <- phyperOverlapSets(gaf, test_genes, background_n)

  shared_desc <- gene_description %>% dplyr::rename(SharedGene = "Gene")
  go_ol$intersecting.genes <- go_ol$intersecting.genes %>%
    left_join(shared_desc, by = "SharedGene")

  if (!is.null(overlaps_file))
    write.csv(go_ol$intersecting.genes, file = overlaps_file, row.names = FALSE)

  scores <- go_ol$pval.scores %>%
    dplyr::rename(GOID           = "Group1",
                  !!group2_label  := "Group2",
                  nGOTermGenes   = "nGroup1",
                  !!ngroup2_label := "nGroup2") %>%
    dplyr::group_by(!!sym(group2_label)) %>%
    dplyr::mutate(FDR = p.adjust(Pval, "BH")) %>%
    dplyr::select(-Enrichment)

  scores <- left_join(go_names, scores, by = "GOID") %>%
    filter(!is.na(nGOTermGenes))

  if (!is.null(pval_filter))
    scores <- scores %>% filter(Pval < pval_filter)
  if (require_shared)
    scores <- scores %>% filter(nSharedGenes > 0)

  if (!is.null(pvals_file))
    write.csv(scores, file = pvals_file, row.names = FALSE)

  return(scores)
}


# run_pairwise_de ---------------------------------------------------------
#' @name run_pairwise_de
#' @description Subsets cells by cluster label, runs findMarkers, writes CSVs.
#' @param combined_sce uncorrected SCE with counts (all sex cells)
#' @param merged_sce batch-corrected SCE supplying $broad and $label colData
#' @param cluster_labels integer vector of label values to keep
#' @param gene_description gene descriptions df
#' @param results_dir output directory
#' @param file_prefix prefix for output file names
#' @param group_col colData column used as grouping variable (default "broad")
#' @param block_col colData column used as blocking variable (default "sample")
#' @return list with $all (all FDR<0.01 genes) and $top100

run_pairwise_de <- function(combined_sce, merged_sce, cluster_labels,
                             gene_description, results_dir, file_prefix,
                             group_col = "broad", block_col = "sample") {
  require(scran)
  require(dplyr)

  colData(combined_sce)$broad <- colData(merged_sce)$broad
  colData(combined_sce)$label <- colData(merged_sce)$label

  keep <- as.data.frame(colData(combined_sce)) %>%
    filter(label %in% cluster_labels) %>%
    pull(barcode)

  sce_sub <- combined_sce[, keep]
  colData(sce_sub) <- droplevels(colData(sce_sub))

  scatter_markers <- findMarkers(sce_sub,
                                  sce_sub[[group_col]],
                                  test.type  = "t",
                                  block      = sce_sub[[block_col]],
                                  direction  = "up",
                                  pval.type  = "some")

  saveRDS(scatter_markers,
          file = paste0(results_dir, "/", file_prefix, "_DEgenes.Rds"))

  scatter_summary <- list()
  for (n in names(scatter_markers)) {
    scatter_markers[[n]]$cluster <- n
    scatter_markers[[n]]$Gene    <- rownames(scatter_markers[[n]])
    scatter_summary[[n]] <- as.data.frame(scatter_markers[[n]]) %>%
      select(cluster, Gene, p.value, FDR, summary.logFC)
  }

  markers_df <- do.call(rbind, scatter_summary) %>%
    left_join(gene_description, by = "Gene") %>%
    filter(FDR < .01)

  write.csv(markers_df,
            file = paste0(results_dir, "/", file_prefix, "_DEgenes.csv"))

  top100 <- markers_df %>%
    group_by(cluster) %>%
    slice_min(FDR, n = 100, with_ties = TRUE) %>%
    mutate(p.value       = formatC(signif(p.value, digits = 3)),
           FDR           = formatC(signif(FDR, digits = 3)),
           summary.logFC = round(summary.logFC, digits = 3))

  return(list(all = markers_df, top100 = top100))
}


# plot_marker_overlap_heatmap ---------------------------------------------
#' @name plot_marker_overlap_heatmap
#' @description Hypergeometric overlap between known cell-type markers and
#'   cluster enriched genes; draws and optionally saves a heatmap.
#' @param known_markers two-column df: Gene, CellType
#' @param test_genes two-column df: Gene, cluster (or other group label)
#' @param background_n number of genes in the universe
#' @param gene_description gene descriptions df
#' @param overlaps_file path for intersecting genes CSV; NULL to skip
#' @param pvals_file path for overlap p-values CSV; NULL to skip
#' @param plot_file path to save the PNG; NULL to skip
#' @return ggplot object

plot_marker_overlap_heatmap <- function(known_markers, test_genes,
                                         background_n, gene_description,
                                         overlaps_file = NULL,
                                         pvals_file    = NULL,
                                         plot_file     = NULL) {
  require(dplyr)
  require(ggplot2)
  require(tidyr)

  cluster_ol <- phyperOverlapSets(known_markers, test_genes, background_n)

  shared_desc <- gene_description %>% dplyr::rename(SharedGene = "Gene")
  cluster_ol$intersecting.genes <- cluster_ol$intersecting.genes %>%
    left_join(shared_desc, by = "SharedGene")

  if (!is.null(overlaps_file))
    write.csv(cluster_ol$intersecting.genes, file = overlaps_file, row.names = FALSE)

  scores <- cluster_ol$pval.scores %>%
    dplyr::rename(CellTypes = "Group1", Clusters = "Group2") %>%
    group_by(CellTypes) %>%
    mutate(FDR = p.adjust(Pval, "BH"))

  if (!is.null(pvals_file))
    write.csv(scores, file = pvals_file, row.names = FALSE)

  FDR_wide <- scores %>%
    mutate(logFDR = -log10(FDR)) %>%
    tidyr::pivot_wider(id_cols = CellTypes, names_from = Clusters, values_from = logFDR)

  CellTypes        <- FDR_wide$CellTypes
  FDR_wide$CellTypes <- NULL
  FDR_mat          <- as.matrix(FDR_wide)
  rownames(FDR_mat) <- CellTypes

  ord_row <- hclust(dist(FDR_mat),    method = "ward.D")$order
  ord_col <- hclust(dist(t(FDR_mat)), method = "ward.D")$order

  scores$CellTypes <- factor(scores$CellTypes, levels = rownames(FDR_mat)[ord_row])
  scores$Clusters  <- factor(scores$Clusters,  levels = colnames(FDR_mat)[ord_col])

  p <- ggplot(scores, aes(Clusters, CellTypes)) +
    geom_tile(aes(fill = -log10(FDR)), color = "light gray") +
    guides(fill = guide_colorbar(title = "-log10(FDR)",
                                  barheight = 10, barwidth = 2,
                                  default.unit = "mm")) +
    scale_fill_gradientn(
      colors = c("white", "#fd9668", "#f1605d", "#de4968", "#9e2f7f"),
      guide  = "colorbar",
      limits = c(0, ceiling(-log10(min(scores$FDR))))
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
    theme(text = element_text(size = 7), axis.text = element_text(size = 7)) +
    geom_text(aes(label = nSharedGenes), size = 1) +
    coord_fixed()

  if (!is.null(plot_file))
    ggsave(p, file = plot_file)

  return(p)
}


# plot_violin_goi ---------------------------------------------------------
#' @name plot_violin_goi
#' @description Stacked violin + boxplot for a set of genes across clusters.
#'   Cluster labels are remapped from numeric `label` column via the `anno`
#'   vector. Significance stars are overlaid from a marker table. Genes are
#'   ordered by hierarchical clustering of mean expression across clusters.
#' @param sce SingleCellExperiment with logcounts and a numeric `label` colData column
#' @param markers data frame with columns WangGeneID, cluster (e.g. "cluster_1"), FDR
#' @param goi character vector of WangGeneIDs to plot
#' @param gene_description gene descriptions df with Gene and WangGeneID columns
#' @param anno character vector mapping numeric cluster index to cell type name
#' @param palette character vector of colors (one per cluster)
#' @param label_col colData column for x-axis grouping; accepts a column name
#'   (character) or column number (integer index into colData). Default "label".
#' @param level_order integer vector of cluster numbers in the same order as anno
#'   (e.g. c(1:4, 6, 8:15, 7, 5, 17, 16, 18, 19)). Required when anno is not
#'   indexed directly by cluster number.
#' @param title optional plot title

plot_violin_goi <- function(sce, markers, goi, gene_description, anno, palette,
                             label_col = "label", level_order = NULL,
                             ncol = 1, title = NULL) {

  fdr_to_stars <- function(fdr) {
    dplyr::case_when(
      fdr < 0.05 ~ "*",
      TRUE       ~ ""
    )
  }

  # build a named lookup: cluster number → cell type label
  if (!is.null(level_order)) {
    anno_lookup <- setNames(anno, as.character(level_order))
  } else {
    anno_lookup <- setNames(anno, as.character(seq_along(anno)))
  }

  if (is.numeric(label_col)) label_col <- names(colData(sce))[label_col]

  # append cell counts to labels
  count_tab   <- table(as.character(colData(sce)[[label_col]]))
  anno_lookup <- setNames(
    paste0(anno_lookup, " (n=", count_tab[names(anno_lookup)], ")"),
    names(anno_lookup)
  )
  anno_n <- unname(anno_lookup[as.character(if (!is.null(level_order)) level_order else seq_along(anno))])

  .to_celltype <- function(x) {
    num <- as.character(as.numeric(gsub("cluster_", "", as.character(x))))
    unname(anno_lookup[num])
  }

  goi_loc <- gene_description %>%
    filter(WangGeneID %in% goi) %>%
    pull(Gene)

  plot_df <- plotExpression(
    sce,
    features      = goi,
    scattermore   = TRUE,
    point_size    = 0,
    swap_rownames = "WangGeneID",
    x             = label_col,
    ncol          = 1
  )$data %>%
    mutate(
      X = .to_celltype(X),
      X = factor(X, levels = anno_n)
    )

  y_max_df <- plot_df %>%
    group_by(Feature, X) %>%
    summarise(y_pos = max(Y) * 1.1, .groups = "drop")

  expand_df <- plot_df %>%
    group_by(Feature) %>%
    summarise(y_expand = max(Y) * 1.4, .groups = "drop") %>%
    mutate(X = anno_n[1])

  sig_df <- markers %>%
    filter(WangGeneID %in% goi) %>%
    mutate(
      stars   = fdr_to_stars(FDR),
      Feature = WangGeneID,
      X       = .to_celltype(cluster)
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

  cluster_colors <- setNames(palette[seq_along(anno_n)], anno_n)

  ggplot(plot_df, aes(x = X, y = Y)) +
    geom_violin(aes(fill = X, color = after_scale(fill)),
                linewidth = 0.3, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.12, fill = "white", color = "grey40",
                 linewidth = 0.3, outlier.shape = NA) +
    geom_blank(data = expand_df, aes(x = X, y = y_expand), inherit.aes = FALSE) +
    geom_text(data = sig_df, aes(x = X, y = y_pos, label = stars),
              inherit.aes = FALSE, size = 2.8, vjust = 0, fontface = "bold",
              color = "#2e2e2e") +
    scale_fill_manual(values = cluster_colors, name = "Cell Type") +
    coord_cartesian(clip = "off", ylim = c(0, NA)) +
    facet_wrap(~ Feature, ncol = ncol, scales = "free_y") +
    labs(x = NULL, y = "Log-normalized expression", title = title) +
    theme_classic(base_size = 8) +
    theme(
      strip.background = element_rect(fill = "#dce8ef", color = NA),
      strip.text       = element_text(color = "#2e4a5a", face = "bold", size = 8),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.text.y      = element_text(color = "#3a3a3a", size = 8),
      axis.title       = element_text(color = "#2e2e2e", size = 8),
      axis.line        = element_line(color = "#aaaaaa", linewidth = 0.4),
      axis.ticks       = element_line(color = "#aaaaaa", linewidth = 0.4),
      panel.spacing    = unit(0.8, "lines"),
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 8),
      legend.key.size  = unit(0.4, "cm")
    )
}


# plot_violin_goi_named ---------------------------------------------------
#' @name plot_violin_goi_named
#' @description Stacked violin + boxplot for a set of genes across clusters
#'   where the cluster column already contains cell type names (no numeric
#'   remapping needed). Genes are ordered by hierarchical clustering of mean
#'   expression. Significance stars from a marker table are overlaid.
#' @param sce SingleCellExperiment with logcounts
#' @param markers data frame with columns WangGeneID, a cluster column matching
#'   cluster_col, and FDR
#' @param goi character vector of WangGeneIDs to plot
#' @param gene_description gene descriptions df with Gene and WangGeneID columns
#' @param anno character vector of cell type names (defines x-axis order and colors)
#' @param palette character vector of colors (one per cell type)
#' @param cluster_col colData column in sce holding cell type names (default "label")
#' @param marker_col column in markers df holding cell type names; defaults to
#'   cluster_col when NULL (set explicitly when the two column names differ)
#' @param title optional plot title

plot_violin_goi_named <- function(sce, markers, goi, gene_description, anno,
                                   palette, cluster_col = "label",
                                   marker_col = NULL, ncol = 1, title = NULL) {
  if (is.null(marker_col)) marker_col <- cluster_col

  fdr_to_stars <- function(fdr) {
    dplyr::case_when(
      fdr < 0.05 ~ "*",
      TRUE       ~ ""
    )
  }

  # append cell counts to anno labels
  count_tab <- table(as.character(colData(sce)[[cluster_col]]))
  anno_n    <- paste0(anno, " (n=", count_tab[anno], ")")
  names(anno_n) <- anno

  plot_df <- plotExpression(
    sce,
    features      = goi,
    scattermore   = TRUE,
    point_size    = 0,
    swap_rownames = "WangGeneID",
    x             = cluster_col,
    ncol          = 1
  )$data %>%
    mutate(
      X = anno_n[as.character(X)],
      X = factor(X, levels = anno_n)
    )

  y_max_df <- plot_df %>%
    group_by(Feature, X) %>%
    summarise(y_pos = max(Y) * 1.1, .groups = "drop")

  expand_df <- plot_df %>%
    group_by(Feature) %>%
    summarise(y_expand = max(Y) * 1.4, .groups = "drop") %>%
    mutate(X = anno_n[1])

  sig_df <- markers %>%
    filter(WangGeneID %in% goi) %>%
    mutate(
      stars   = fdr_to_stars(FDR),
      Feature = WangGeneID,
      X       = factor(anno_n[as.character(.data[[marker_col]])], levels = anno_n)
    ) %>%
    filter(stars != "") %>%
    left_join(y_max_df, by = c("Feature", "X"))

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

  cluster_colors <- setNames(palette[seq_along(anno_n)], anno_n)

  ggplot(plot_df, aes(x = X, y = Y)) +
    geom_violin(aes(fill = X, color = after_scale(fill)),
                linewidth = 0.3, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.12, fill = "white", color = "grey40",
                 linewidth = 0.3, outlier.shape = NA) +
    geom_blank(data = expand_df, aes(x = X, y = y_expand), inherit.aes = FALSE) +
    geom_text(data = sig_df, aes(x = X, y = y_pos, label = stars),
              inherit.aes = FALSE, size = 2.8, vjust = 0, fontface = "bold",
              color = "#2e2e2e") +
    scale_fill_manual(values = cluster_colors, name = "Cell Type") +
    coord_cartesian(clip = "off", ylim = c(0, NA)) +
    facet_wrap(~ Feature, ncol = ncol, scales = "free_y") +
    labs(x = NULL, y = "Log-normalized expression", title = title) +
    theme_classic(base_size = 8) +
    theme(
      strip.background = element_rect(fill = "#dce8ef", color = NA),
      strip.text       = element_text(color = "#2e4a5a", face = "bold", size = 8),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.text.y      = element_text(color = "#3a3a3a", size = 8),
      axis.title       = element_text(color = "#2e2e2e", size = 8),
      axis.line        = element_line(color = "#aaaaaa", linewidth = 0.4),
      axis.ticks       = element_line(color = "#aaaaaa", linewidth = 0.4),
      panel.spacing    = unit(0.8, "lines"),
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 8),
      legend.key.size  = unit(0.4, "cm")
    )
}


# plot_violin_goi_single --------------------------------------------------
#' @name plot_violin_goi_single
#' @description Horizontal (flipped) violin + boxplot for a single gene across
#'   clusters. Cluster labels are remapped from numeric `label` via `anno`.
#'   Significance stars from a marker table are overlaid.
#' @param sce SingleCellExperiment with logcounts and a numeric `label` colData column
#' @param markers data frame with columns WangGeneID, cluster (e.g. "cluster_1"), FDR
#' @param goi single WangGeneID string to plot
#' @param gene_description gene descriptions df with Gene and WangGeneID columns
#' @param anno character vector mapping numeric cluster index to cell type name
#' @param palette character vector of colors (one per cluster)
#' @param level_order integer vector of cluster numbers in the same order as anno
#'   (e.g. c(1:4, 6, 8:15, 7, 5, 17, 16, 18, 19)). Required when anno is not
#'   indexed directly by cluster number.
#' @param title optional plot title

plot_violin_goi_single <- function(sce, markers, goi, gene_description, anno,
                                    palette, level_order = NULL, title = NULL) {

  fdr_to_stars <- function(fdr) {
    dplyr::case_when(
      fdr < 0.05 ~ "*",
      TRUE       ~ ""
    )
  }

  if (!is.null(level_order)) {
    anno_lookup <- setNames(anno, as.character(level_order))
  } else {
    anno_lookup <- setNames(anno, as.character(seq_along(anno)))
  }

  count_tab   <- table(as.character(colData(sce)[["label"]]))
  anno_lookup <- setNames(
    paste0(anno_lookup, " (n=", count_tab[names(anno_lookup)], ")"),
    names(anno_lookup)
  )
  anno_n <- unname(anno_lookup[as.character(if (!is.null(level_order)) level_order else seq_along(anno))])

  .to_celltype <- function(x) {
    num <- as.character(as.numeric(gsub("cluster_", "", as.character(x))))
    unname(anno_lookup[num])
  }

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
      X = .to_celltype(X),
      X = factor(X, levels = rev(anno_n))
    )

  y_max_df <- plot_df %>%
    group_by(Feature, X) %>%
    summarise(y_pos = max(Y) * 1.1, .groups = "drop")

  expand_df <- plot_df %>%
    group_by(Feature) %>%
    summarise(y_expand = max(Y) * 1.4, .groups = "drop") %>%
    mutate(X = anno_n[length(anno_n)])

  sig_df <- markers %>%
    filter(WangGeneID %in% goi) %>%
    mutate(
      stars   = fdr_to_stars(FDR),
      Feature = WangGeneID,
      X       = .to_celltype(cluster)
    ) %>%
    filter(stars != "") %>%
    left_join(y_max_df, by = c("Feature", "X"))

  cluster_colors <- setNames(palette[seq_along(anno_n)], anno_n)

  ggplot(plot_df, aes(x = X, y = Y)) +
    geom_violin(aes(fill = X, color = after_scale(fill)),
                linewidth = 0.3, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.12, fill = "white", color = "grey40",
                 linewidth = 0.3, outlier.shape = NA) +
    geom_blank(data = expand_df, aes(x = X, y = y_expand), inherit.aes = FALSE) +
    geom_text(data = sig_df, aes(x = X, y = y_pos, label = stars),
              inherit.aes = FALSE, size = 2.8, vjust = 0, fontface = "bold",
              color = "#2e2e2e") +
    scale_fill_manual(values = cluster_colors, name = "Cell Type") +
    coord_flip(clip = "off", ylim = c(0, NA)) +
    labs(x = NULL, y = "Log-normalized expression", title = title) +
    theme_classic(base_size = 8) +
    theme(
      axis.text.x     = element_text(color = "#3a3a3a", size = 8),
      axis.ticks.x    = element_blank(),
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.title      = element_text(color = "#2e2e2e", size = 8),
      axis.line       = element_line(color = "#aaaaaa", linewidth = 0.4),
      axis.ticks      = element_line(color = "#aaaaaa", linewidth = 0.4),
      legend.title    = element_text(size = 8),
      legend.text     = element_text(size = 8),
      legend.key.size = unit(0.4, "cm")
    ) +
    guides(fill = guide_legend(reverse = TRUE))
}


# plot_expr_proportion ----------------------------------------------------
#' @name plot_expr_proportion
#' @description Bar chart showing the proportion of cells expressing a gene
#'   (above a threshold) per cell type, faceted by batch/sample. Useful for
#'   checking whether expression is sample-specific.
#' @param sce SingleCellExperiment object
#' @param gene rowname in sce to plot
#' @param celltype_col colData column with cell type labels (default "celltype")
#' @param batch_col colData column with batch/sample labels (default "batch")
#' @param assay_name assay to use for expression values (default "counts")
#' @param threshold expression threshold above which a cell is "expressing"
#'   (default 0)
#' @param ncol number of columns in the facet grid (default 3)

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
    geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)),
              vjust = -0.4, size = 3) +
    facet_wrap(~ celltype, ncol = ncol) +
    scale_y_continuous(labels = scales::percent_format(),
                       limits = c(0, 1.20),
                       breaks = seq(0, 1, by = 0.25)) +
    labs(title = paste("Proportion of cells expressing", gene),
         x     = "Batch",
         y     = "% cells expressing",
         fill  = "Batch") +
    theme_bw() +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1),
          legend.position = "none")
}


# plot_pseudotime_paths ---------------------------------------------------
#' @name plot_pseudotime_paths
#' @description GAM-smoothed expression curves along pseudotime, one line per
#'   trajectory path, faceted by gene. FDR values from a pseudotime DE results
#'   table are annotated in the top-right of each facet.
#' @param sce SingleCellExperiment with logcounts
#' @param pseudo_mat matrix from pathStat(); rows = cells, columns = path endpoints
#' @param goi character vector of WangGeneIDs to plot
#' @param gene_description gene descriptions df with Gene and WangGeneID columns
#' @param tradeseq_res data frame of pseudotime DE results with columns
#'   WangGeneID, Path, FDR (e.g. from testPseudotime + formatting)
#' @param paths character vector of path (column) names from pseudo_mat to include
#'   (default c("Oocytes", "SPC.cytokinetic"))
#' @param title optional plot title

plot_pseudotime_paths <- function(sce, pseudo_mat, goi, gene_description,
                                   tradeseq_res,
                                   paths = c("Oocytes", "SPC.cytokinetic"),
                                   title = NULL) {

  goi_loc <- gene_description %>%
    filter(WangGeneID %in% goi) %>%
    pull(Gene)

  expr_mat <- as.matrix(logcounts(sce)[goi_loc, , drop = FALSE])
  rownames(expr_mat) <- gene_description %>%
    filter(Gene %in% goi_loc) %>%
    select(Gene, WangGeneID) %>%
    deframe() %>%
    { .[goi_loc] }

  pseudo_long <- pseudo_mat[, paths, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("Cell") %>%
    tidyr::pivot_longer(-Cell, names_to = "Path", values_to = "Pseudotime") %>%
    filter(!is.na(Pseudotime))

  expr_long <- as.data.frame(t(expr_mat)) %>%
    rownames_to_column("Cell") %>%
    tidyr::pivot_longer(-Cell, names_to = "Gene", values_to = "Expression")

  plot_df <- pseudo_long %>%
    left_join(expr_long, by = "Cell") %>%
    filter(!is.na(Expression))

  fdr_labels <- tradeseq_res %>%
    filter(WangGeneID %in% goi, Path %in% paths) %>%
    mutate(label = paste0(Path, ": ", formatC(FDR, format = "e", digits = 2))) %>%
    group_by(WangGeneID) %>%
    summarise(fdr_text = paste(label, collapse = "\n"), .groups = "drop") %>%
    dplyr::rename(Gene = WangGeneID)

  path_colors <- setNames(c("#E377C2", "#1F77B4"), paths)

  ggplot(plot_df, aes(x = Pseudotime, y = Expression, color = Path, fill = Path)) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                se = TRUE, alpha = 0.15, linewidth = 0.8) +
    geom_text(data = fdr_labels,
              aes(x = Inf, y = Inf, label = fdr_text),
              inherit.aes = FALSE,
              hjust = 1.05, vjust = 1.2, size = 2.5, color = "grey30",
              lineheight = 1.4) +
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


# plot_tsne ---------------------------------------------------------------
#' @name plot_tsne
#' @description ggplot2 scatter plot of TSNE coordinates with optional cluster
#'   text labels. Mirrors plotTSNE() aesthetics (black background, no axes).
#' @param df data frame with columns tsne1, tsne2, and any metadata columns
#' @param color_col column name in df to map to point color
#' @param label_col column name in df with cell type labels to display at
#'   cluster medians. If NULL, no text is drawn.
#' @param palette optional named character vector of colors passed to
#'   scale_color_manual(). If NULL, ggplot2 default colors are used.
#' @param point_size point size (default 0.05)
#' @param text_color color for cluster label text (default "red")
#' @param text_size size for cluster label text (default 3)
#' @param bg_color background color for panel, plot, and legend (default "black");
#'   legend and title text automatically use black or white based on luminance
#' @param color_levels optional character vector of levels for color_col, controls
#'   legend order (e.g. anno[level_order])
#' @param title optional plot title

plot_tsne <- function(df,
                      color_col,
                      label_col    = NULL,
                      palette      = NULL,
                      point_size   = 0.05,
                      text_color   = "red",
                      text_size    = 3,
                      bg_color     = "black",
                      color_levels = NULL,
                      title        = NULL) {

  # pick legible text color based on bg luminance (WCAG relative luminance)
  .contrast_color <- function(hex) {
    rgb  <- col2rgb(hex) / 255
    lin  <- ifelse(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
    lum  <- 0.2126 * lin[1] + 0.7152 * lin[2] + 0.0722 * lin[3]
    if (lum > 0.179) "black" else "white"
  }
  fg_color <- .contrast_color(bg_color)

  if (!is.null(color_levels))
    df[[color_col]] <- factor(df[[color_col]], levels = color_levels)

  p <- ggplot(df, aes(x = TSNE1, y = TSNE2, color = .data[[color_col]])) +
    geom_point(size = point_size, stroke = 0) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    labs(color = color_col, title = title) +
    theme_classic() +
    theme(
      panel.background  = element_rect(fill = bg_color, color = NA),
      plot.background   = element_rect(fill = bg_color, color = NA),
      axis.title        = element_blank(),
      axis.text         = element_blank(),
      axis.ticks        = element_blank(),
      axis.line         = element_blank(),
      legend.background = element_rect(fill = bg_color),
      legend.text       = element_text(color = fg_color),
      legend.title      = element_text(color = fg_color),
      plot.title        = element_text(color = fg_color)
    )

  if (!is.null(palette)) {
    if (!is.null(color_levels) && is.null(names(palette)))
      palette <- setNames(palette[seq_along(color_levels)], color_levels)
    p <- p + scale_color_manual(values = palette,
                                breaks = if (!is.null(color_levels)) color_levels else waiver())
  }

  if (!is.null(label_col)) {
    label_df <- df %>%
      group_by(.data[[label_col]]) %>%
      summarise(TSNE1 = median(TSNE1), TSNE2 = median(TSNE2), .groups = "drop") %>%
      rename(.label = all_of(label_col))

    p <- p +
      geom_point(
        data        = label_df,
        aes(x = TSNE1, y = TSNE2),
        color       = "white",
        size        = text_size * 2,
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      geom_text(
        data        = label_df,
        aes(x = TSNE1, y = TSNE2, label = .label),
        color       = text_color,
        size        = text_size,
        fontface    = "bold",
        inherit.aes = FALSE
      )
  }

  p
}
