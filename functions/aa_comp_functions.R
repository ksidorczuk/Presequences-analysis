#' Get heatmap of amino acid composition of datasets
#' 
#' This function uses a data frame of amino acid composition to plot
#' a heatmap with a dendrogram.
#' @param aa_comp_all amino acid composition of whole datasets
get_aa_comp_heatmap <- function(aa_comp_all) {
  aa_comp_heatmap_dat <- aa_comp_all %>% 
    pivot_wider(names_from = "aa", values_from = "Freq", values_fill = 0) 
  
  clustering_datasets <- as.dendrogram(hclust(dist(as.matrix(aa_comp_heatmap_dat[,2:21]))))
  dataset_order <- order.dendrogram(clustering_datasets)
  
  clustering_aa <- as.dendrogram(hclust(dist(t(as.matrix(aa_comp_heatmap_dat[,2:21])))))
  aa_order <- order.dendrogram(clustering_aa)
  
  dendro <- clustering_datasets %>%
    dendro_data %>%
    segment %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_segment() +
    scale_y_continuous("") +
    scale_x_discrete("",
                     limits = factor(1L:nobs(clustering_datasets))) + 
    theme_void() + 
    coord_flip()
  
  aa_comp_clustered_dat <- pivot_longer(aa_comp_heatmap_dat, 2:21, names_to = "Amino acid", values_to = "Frequency")
  aa_comp_clustered_dat[["dataset"]] <- factor(aa_comp_clustered_dat[["dataset"]],
                                               levels = aa_comp_heatmap_dat[["dataset"]][dataset_order], ordered = TRUE)
  aa_comp_clustered_dat[["Amino acid"]] <- factor(aa_comp_clustered_dat[["Amino acid"]],
                                                  levels = unique(aa_comp_clustered_dat[["Amino acid"]])[aa_order], ordered = TRUE)
  
  heatmap <- aa_comp_clustered_dat %>% 
    ggplot(aes(x = `Amino acid`, y = dataset)) +
    geom_tile(aes(fill = Frequency)) +
    scale_fill_gradientn(colors = c("#ffffff", "#ffe96b", "#ff4242", "#630000"), values = rescale(c(0, 0.03, 0.08, 0.14), to = c(0, 1))) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.key.width = unit(2, "lines")) +
    ylab("Data set")
  
  max_height <- unit.pmax(ggplotGrob(heatmap)[["heights"]],
                          ggplotGrob(dendro)[["heights"]])
  
  grob_list <- list(heatmap = ggplotGrob(heatmap),
                    dendrogram = ggplotGrob(dendro))
  grob_list[["heatmap"]][["heights"]] <-
    grob_list[["dendrogram"]][["heights"]] <-
    max_height
  
  grid.arrange(grob_list[["heatmap"]], grob_list[["dendrogram"]], widths = c(0.8, 0.2))
}


#' Test differences in amino acid composition
#' 
#' This function performs Mann-Whitney test to compare amino acid composition
#' of peptides from different datasets.
#' @param aa_comp_peptides amino acid composition of peptides from each dataset
#' obtained using \code{\link{calculate_aa_comp_peptides}}
#' @param combns a list of dataset pairs that will be used for testing
#' @return a data frame containing names of the compared datasets, amino acid, 
#' p-value adjusted for multiple comparisons and information if the difference
#' is statistically significant (p < 0.05).
get_statistical_analysis_aa_comp <- function(aa_comp_peptides, combns) {

  lapply(seq_along(combns), function(ith_combn) {
    test_dat <- filter(aa_comp_peptides, dataset  %in% combns[[ith_combn]])
    lapply(unique(aa_comp_peptides[["Amino acid"]]), function(ith_aa) {
      data.frame(dataset1 = combns[[ith_combn]][1], 
                 dataset2 = combns[[ith_combn]][2],
                 aa = ith_aa,
                 pval = wilcox.test(x = filter(test_dat, `Amino acid` == ith_aa, dataset == combns[[ith_combn]][1])[["Frequency"]],
                                    y = filter(test_dat, `Amino acid` == ith_aa, dataset == combns[[ith_combn]][2])[["Frequency"]],
                                    exact = FALSE)[["p.value"]])
    }) %>% bind_rows() 
  }) %>% bind_rows() %>%
    mutate(pval_adjusted = p.adjust(pval)) %>%
    select(-pval) %>%
    group_by(aa, dataset1, dataset2) %>%
    dplyr::summarise(is_signif = pval_adjusted < 0.05) 
}
