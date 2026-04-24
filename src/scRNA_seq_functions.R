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

# modify_vlnplot ----------------------------------------------------------
#' @name modify_vlnplot
#' @description from https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
#' @description remove the x-axis text and tick marks from a Seurat Vln plot
#' @param obj seurat object
#' @param feature  one genes to be plotted
#' @param pt.size size of points on Vln plot
#' @param plot.margin to adjust the white space between each plot
#' @param pass any arguments to VlnPlot in Seurat

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  require(Seurat)
  require(patchwork)
  
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab("norm expr") + ggtitle(feature) + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = 8), 
          axis.text.y = element_text(size = 8),
          plot.title = element_text(size = 8),
          plot.margin = plot.margin
          ) +
    ## add box plot (steph edit)
    geom_boxplot(width=0.1, color="black") +
    theme(legend.position = 'none')
  return(p)
}

# extract_max_y -------------------------------------------------------------
#' @name extract_max_y
#' @description given a ggplot, extract the max value of the y axis
#' @param p a ggplot object

extract_max_y<- function(p){
  require(ggplot2)
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

# StackedVlnPlot_samey ----------------------------------------------------
#' @name StackedVlnPlot_samey
#' @description make a stacked vln plot with the same y-axis for each gene
#' @param obj seurat object
#' @param features  charcter vector of genes to be plotted
#' @param pt.size size of points on Vln plot
#' @param plot.margin to adjust the white space between each plot
#' @param pass any arguments to VlnPlot in Seurat
 
StackedVlnPlot_samey<- function(obj, 
                                features,
                                pt.size = 0, 
                                plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                                ...) {
  
  require(Seurat)
  require(patchwork)
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 30,  hjust = 1, vjust = 1), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max_y)
  ## fill ymaxs with the max of all of the plots (steph edit)
  same_y = rep(max(ymaxs), length(ymaxs))
  plot_list<- purrr::map2(plot_list, same_y, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

# StackedVlnPlot_samey ----------------------------------------------------
#' @name StackedVlnPlot_samey
#' @description make a stacked vln plot with the max-y for each gene
#' @param obj seurat object
#' @param features  charcter vector of genes to be plotted
#' @param pt.size size of points on Vln plot
#' @param plot.margin to adjust the white space between each plot
#' @param pass any arguments to VlnPlot in Seurat

StackedVlnPlot<- function(obj, features, plot_title,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  require(Seurat)
  require(patchwork)
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 30, hjust = 1, vjust = 1, size = 8), axis.ticks.x = element_line()) 
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max_y)
  same_y = rep(max(ymaxs), length(ymaxs))
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  p = p + plot_annotation(title = plot_title) &
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(face="bold")) + 
    theme(text = element_text(size = 8))
  return(p)
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


# plot_trajectory_genes ---------------------------------------------------
# line plot of average gene expression across a trajectory
plot_trajectory_genes <- function(sce, genes, cluster_order, 
                                  cluster_col, title, assay = "logcounts") {
  
  avg_exp <- get_avg_expression(sce, genes, cluster_col, assay)
  
  wang_genes = 
    gene_description %>%
    filter(Gene %in% genes) %>%
    pull(WangGeneID)
  
  # Reorder to trajectory order
  valid_clusters <- cluster_order[cluster_order %in% colnames(avg_exp)]
  avg_exp <- avg_exp[, valid_clusters, drop = FALSE]
  
  # Scale each gene 0-1
  avg_exp_scaled <- t(apply(avg_exp, 1, function(x) {
    rng <- max(x) - min(x)
    if (rng == 0) return(rep(0, length(x)))
    (x - min(x)) / rng
  }))
  
  # Long format for ggplot
  df <- as.data.frame(avg_exp_scaled)
  df$gene <- wang_genes
  df_long <- pivot_longer(df, -gene, 
                          names_to = "cluster", 
                          values_to = "expression")
  df_long$cluster <- factor(df_long$cluster, levels = valid_clusters)
  
  library(Polychrome)
  n_clusters <- length(wang_genes)
  color_pal <- as.vector(kelly.colors(n_clusters + 1))[-1]  # up to 22 colors
  
  
  ggplot(df_long, aes(x = cluster, y = expression, 
                      color = gene, group = gene)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    labs(title = title,
         x = "Cluster (ordered by trajectory)",
         y = "Scaled mean expression (0-1)") +
    scale_color_manual(values = color_pal)
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
