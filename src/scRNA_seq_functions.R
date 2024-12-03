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