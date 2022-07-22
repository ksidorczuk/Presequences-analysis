#' Calculate ngram frequency
#' 
#' Calculates frequency of ngrams in a list of sequence
#' datasets. Motifs absent in a given dataset will have
#' frequency 0. Frequency is calculated for each dataset 
#' as a fraction of sequences in which a motif is present.
#' 
#' @param datasets a named list of sequence datasets 
#' @param k a vector of k values (ngram size)
#' @param kmer_gaps a list of ngram gap sizes 
#' 
#' @return a data frame with datasets as rows and ngrams
#' as columns with additional column specifying dataset name
calculate_ngram_freq <- function(datasets, k, kmer_gaps, alphabet = NULL) {
  freqs <- lapply(names(datasets), function(ith_dataset) {
    x <- calculate_binary_ngram_matrix(datasets[[ith_dataset]],
                                       k = k,
                                       kmer_gaps_list = kmer_gaps) 
    if(!is.null(alphabet)) {
      x <- degenerate_ngrams(x, alphabet, binarize = TRUE)
    }
    freq <- colSums(x)/length(datasets[[ith_dataset]])
    t(freq) %>% 
      as.data.frame() %>% 
      mutate(dataset = ith_dataset)
  }) %>% bind_rows() 
  freqs[is.na(freqs)] <- 0
  freqs
}

#' Locate motifs
#' 
#' Identifies motif start and end positions within sequences.
#' 
#' @param dataset a list of sequences 
#' @param motifs a vector of motifs to find (written as
#' regular expressions)
#' 
#' @return a data frame with the following columns: start/end 
#' positions of a motif, dataset name, sequence name and motif
locate_motifs <- function(dataset, motifs, alphabet = NULL) {
  pblapply(motifs, cl = 4, function(ith_motif) {
    if(is.null(alphabet)) {
      ds <- unlist(unname(dataset), recursive = FALSE)
    } else {
      ds <- sapply(unlist(unname(dataset), recursive = FALSE), function(i) degenerate(i, alphabet))
    }
    lapply(names(ds), function(ith_seq) {
      str_locate_all(paste0(ds[[ith_seq]], collapse = ""), ith_motif) %>%
        as.data.frame() %>%
        mutate(dataset = names(dataset),
               seq = ith_seq,
               motif = ith_motif)
    }) %>% bind_rows()
  }) %>% bind_rows()
}


#' Get motif plot data
#' 
#' Get long motif data for plots
#' 
#' @param motif_freqs a data frame with motif frequencies in
#' datasets (obtained with \code{calculate_ngram_freq})
#' @param gapped a logical indicating if motifs are gapped
#' (needed for motif decoding into regular expression)
#' 
#' @return a long data frame with the following columns: dataset
#' name, motif, frequency of the motif
get_motif_plot_dat <- function(motif_freqs) {
  colnames(motif_freqs) <- sapply(colnames(motif_freqs), function(i) {
    ifelse(i == "dataset", "dataset", gsub("_", ".", decode_ngrams(i), fixed = TRUE))
  })
  pivot_longer(motif_freqs, colnames(motif_freqs)[which(colnames(motif_freqs) != "dataset")], 
               names_to = "motif", values_to = "frequency")
}


#' Get most frequent motifs for a dataset
#' 
#' Extracts a given number of most frequent
#' motifs in a given dataset.
#' 
#' @param motif_plot_data motif data for plots,
#' obtained with \code{get_motif_plot_dat}
#' @param dataset_name name of a dataset for which
#' motifs are to be extracted
#' @param n_motifs number of motifs to extract
#' 
#' @return a vector of the most frequent motifs
get_dataset_most_frequent_motifs <- function(motif_plot_data, dataset_name, n_motifs) {
  motif_plot_data %>% 
    filter(dataset == dataset_name) %>%  
    arrange(desc(frequency)) %>% 
    head(n = n_motifs) %>% 
    .[["motif"]]
}


#' Plot most frequent motifs in a dataset
#' 
#' Plots a bar plot with frequencies of selected motifs
#' in all datasets. Plot is faceted using ordered motifs,
#' so they are ordered according to the frequency in the
#' one dataset.
#' 
#' @param motif_plot_data
#' @param selected_motifs a vector of motifs to plot 
#' (their order will be used to order facets)
#' @param colors a named vector of colors for each dataset
#' @param title plot title, e.g. specifying dataset for
#' which the most frequent motifs were used
#' 
#' @return a faceted bar plot with motif frequencies 
plot_most_frequent_motifs <- function(motif_plot_data, selected_motifs, colors, title) {
  motif_plot_data %>% 
    filter(motif %in% selected_motifs) %>% 
    mutate(motif = factor(motif, levels = selected_motifs)) %>% 
    ggplot(aes(x = dataset, y = frequency, fill = dataset)) +
    geom_col() +
    facet_wrap(~motif) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_fill_manual("Data set", values = colors) +
    ggtitle(title)
}


#' Plot motif positions within sequences
#' 
#' Plots sequences and motifs within them indicating numbers
#' of motifs in a given position.
#' 
#' @param dataset_list a named list of sequence datasets 
#' @param dataset_name name of a dataset with sequences to plot
#' @param found_motifs a data frame with motif positions, obtained
#' with \code{locate_motifs}
#' 
#' @return a heatmap showing sequences with their lengths and motif positions,
#' faceted by motifs
plot_motif_positions <- function(dataset_list, dataset_name, found_motifs) {
  lens <- data.frame(seq = names(dataset_list[[dataset_name]]),
                     len = lengths(dataset_list[[dataset_name]]))
  
  lens_dat <- left_join(filter(found_motifs, dataset == dataset_name), lens, by = "seq")
  
  plot_dat <- pblapply(1:nrow(lens_dat), cl = 4, function(i) {
    data.frame(dataset = lens_dat[["dataset"]][i],
               seq_name = lens_dat[["seq"]][i],
               motif_pos = lens_dat[["start"]][i]:lens_dat[["end"]][i],
               len = lens_dat[["len"]][i],
               motif = lens_dat[["motif"]][[i]])
  }) %>% bind_rows()
  
  x_len <- lapply(1:nrow(lens), function(i) {
    data.frame(seq = lens[["seq"]][i],
               pos = 1:lens[["len"]][i],
               n = "1")
  }) %>% bind_rows() 
  
  n_ngrams_plot_dat <- plot_dat %>% 
    group_by(dataset, seq_name, len, motif, motif_pos) %>% 
    summarise(n = n()) 
  
  ggplot() +
    geom_tile(data = x_len, aes(x = pos, y = seq, fill = n)) +
    scale_fill_manual(values = "grey70") +
    new_scale_fill() +
    geom_tile(data = n_ngrams_plot_dat, aes(x = motif_pos, y = seq_name, fill = n)) +
    scale_fill_gradient(low = "yellow", high = "red") +
    theme_bw(base_size = 5) +
    facet_wrap(~motif, ncol = 8) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle(dataset_name)
}

#' Plot motif position density
#' 
#' Plots density of sequence positions and motif positions
#' within a given dataset.
#' 
#' @param dataset_list a named list of sequence datasets 
#' @param dataset_name name of a dataset with sequences to plot
#' @param found_motifs a data frame with motif positions, obtained
#' @param density_col color for the violin plots of motif density
#' 
#' @return a density plot of sequence and motif positions, faceted
#' by motifs
plot_motif_position_density <- function(dataset_list, dataset_name, found_motifs, density_col) {
  lens <- data.frame(seq = names(dataset_list[[dataset_name]]),
                     len = lengths(dataset_list[[dataset_name]]))
  
  lens_dat <- left_join(found_motifs, lens, by = "seq")
  
  plot_dat <- pblapply(1:nrow(lens_dat), cl = 4, function(i) {
    data.frame(dataset = lens_dat[["dataset"]][i],
               seq_name = lens_dat[["seq"]][i],
               motif_pos = lens_dat[["start"]][i]:lens_dat[["end"]][i],
               len = lens_dat[["len"]][i],
               motif = lens_dat[["motif"]][[i]])
  }) %>% bind_rows()
  
  n_ngrams_plot_dat <- plot_dat %>% 
    group_by(dataset, seq_name, len, motif, motif_pos) %>% 
    summarise(n = n()) 
  
  pblapply(1:nrow(lens_dat), cl = 4, function(i) {
    data.frame(dataset = lens_dat[["dataset"]][i],
               seq_name = lens_dat[["seq"]][i],
               pos = 1:lens_dat[["len"]][i])
  }) %>% 
    bind_rows() %>% 
    ggplot(aes(x = dataset, y = pos)) +
    geom_violin(fill = "grey90", color = "grey90") +
    coord_flip() +
    geom_violin(data = n_ngrams_plot_dat, aes(x = dataset, y = motif_pos), 
                color = density_col, fill = density_col, alpha = 0.25) +
    facet_wrap(~motif) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle(dataset_name)
}


plot_motif_positions_scaled <- function(dataset_list, dataset_name, found_motifs, type = "point", density_col = "white") {
  lens <- data.frame(seq = names(dataset_list[[dataset_name]]),
                     len = lengths(dataset_list[[dataset_name]]))
  
  lens_dat <- left_join(filter(found_motifs, dataset == dataset_name), lens, by = "seq") %>% 
    mutate(start = start*100/len) 
  
  if(type == "point") {
    ggplot(lens_dat, aes(x = start, y = dataset)) +
      geom_point(shape = 15) +
      facet_wrap(~motif) +
      theme_bw()  +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      ggtitle(dataset_name)
  } else {
    ggplot(lens_dat, aes(x = start, y = dataset)) +
      geom_violin(fill = density_col, alpha = 0.25) +
      facet_wrap(~motif) +
      theme_bw()  +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      ggtitle(dataset_name)
  }
}


calculate_ngram_fractions <- function(datasets, k, kmer_gaps) {
  fracs <- lapply(names(datasets), function(ith_dataset) {
    x <- count_multimers(datasets[[ith_dataset]],
                         k_vector = k,
                         kmer_gaps_list = kmer_gaps, 
                         kmer_alphabet = toupper(colnames(biogram::aaprop)),
                         with_kmer_counts = TRUE) 
    frac <- colSums(as.matrix(x/sum(x)))
    frac %>%
      t() %>% 
      as.data.frame() %>% 
      mutate(dataset = ith_dataset)
  }) %>% bind_rows()
  fracs[is.na(fracs)] <- 0
  fracs
}

select_quipt_informative_motifs <- function(binary_ngrams, dataset1_name, dataset2_name, pval = 0.05) {
  x <- binary_ngrams %>% 
    filter(dataset %in% c(dataset1_name,
                          dataset2_name)) %>% 
    mutate(dataset = ifelse(dataset == dataset1_name, 1, 0)) 
  tf <- test_features(x[["dataset"]], select(x, -dataset))
  filter(as.data.frame(tf), p.value <= pval) %>% 
    mutate(ngram = sapply(.[["ngram"]], function(i) gsub("_", ".", decode_ngrams(i), fixed = TRUE))) %>% 
    arrange(p.value) %>% 
    setNames(c("motif", "p-value", dataset1_name, paste0(dataset2_name, collapse = "/")))
}

select_fcbf_informative_motifs <- function(binary_ngrams, dataset1_name, dataset2_name, min_su = 0.05) {
  x <- binary_ngrams %>% 
    filter(dataset %in% c(dataset1_name,
                          dataset2_name)) %>% 
    mutate(dataset = ifelse(dataset == dataset1_name, 1, 0))
  tar <- x[["dataset"]]
  x <- select(x, -dataset)
  x[sapply(x, is.numeric)] <- lapply(x[sapply(x, is.numeric)], 
                                     as.factor)
  fcbf_res <- tryCatch({
    fcbf(x, as.factor(tar), samples_in_rows = TRUE, minimum_su = min_su)
  },
  error = function(e) NULL)
  if(!is.null(fcbf_res)) {
    if(nrow(fcbf_res) == 1) {
      rownames(fcbf_res) <- colnames(binary_ngrams)[fcbf_res[["index"]][1]]
      fcbf_res
    } else {
      as.data.frame(fcbf_res) %>%
        mutate(motif = sapply(rownames(.), function(i) gsub("_", ".", decode_ngrams(i), fixed = TRUE)))
    }
  } else {
    NULL
  }
}


plot_motif_venn_diagram <- function(datasets, colors, amp = FALSE) {
  ds_list <- if(amp == TRUE) {
    list("cTP-mTP" = filter(datasets, 
                            dataset == "cTP-mTP experimentally verified presequence" & frequency > 0)[["motif"]],
         "cTP" = filter(datasets, 
                        dataset == "cTP experimentally verified presequence" & frequency > 0)[["motif"]],
         "mTP" = filter(datasets, 
                        dataset == "mTP experimentally verified presequence" & frequency > 0)[["motif"]],
         "SP" = filter(datasets, 
                       dataset == "SP experimentally verified presequence" & frequency > 0)[["motif"]],
         "AMP" = filter(datasets, 
                        dataset == "DBAASP AMP max 100 aa" & frequency > 0)[["motif"]])
  } else {
    list("cTP-mTP" = filter(datasets, 
                            dataset == "cTP-mTP experimentally verified presequence" & frequency > 0)[["motif"]],
         "cTP" = filter(datasets, 
                        dataset == "cTP experimentally verified presequence" & frequency > 0)[["motif"]],
         "mTP" = filter(datasets, 
                        dataset == "mTP experimentally verified presequence" & frequency > 0)[["motif"]],
         "SP" = filter(datasets, 
                       dataset == "SP experimentally verified presequence" & frequency > 0)[["motif"]])
  }
  ggVennDiagram(ds_list,
                label_alpha = 0) +
    scale_color_manual("Dataset", values = unname(colors)) +
    scale_fill_continuous(low = "white", high = "tan1")
}


do_statistical_analysis <- function(ngram_binary_matrix) {
  combns <- unique(ngram_binary_matrix[["dataset"]]) %>% 
    combn(., 2, simplify = FALSE)
  colnames(ngram_binary_matrix) <- sapply(colnames(ngram_binary_matrix), function(i) {
    ifelse(i == "dataset", "dataset", gsub("_", ".", decode_ngrams(i), fixed = TRUE))
  })
  lapply(seq_along(combns), function(ith_combn) {
    test_dat <- filter(ngram_binary_matrix, dataset %in% combns[[ith_combn]])
    pblapply(colnames(ngram_binary_matrix)[which(colnames(ngram_binary_matrix) != "dataset")], 
             cl = 3, 
             function(ith_motif) {
               data.frame(comparison = paste0(strsplit(combns[[ith_combn]][1], " ")[[1]][1], "_", 
                                              strsplit(combns[[ith_combn]][2], " ")[[1]][1]),
                          motif = ith_motif,
                          pval = wilcox.test(x = filter(test_dat, dataset == combns[[ith_combn]][1])[[ith_motif]],
                                             y = filter(test_dat, dataset == combns[[ith_combn]][2])[[ith_motif]],
                                             exact = FALSE)[["p.value"]]) 
             }) %>% bind_rows() %>% 
      mutate(pval_adjusted = p.adjust(pval))
  }) %>% bind_rows()
}


plot_clustered_statistical_analysis_res <- function(test_res) {
  test_res_plot_dat <- test_res %>% 
    select(-pval)  %>% 
    group_by(motif, comparison) %>% 
    summarise(is_significant = as.logical(pval_adjusted < 0.05)) %>% 
    pivot_wider(motif, names_from = comparison, values_from = is_significant)
  
  clustering_motifs <- as.dendrogram(hclust(dist(as.matrix(test_res_plot_dat[, 2:7]))))
  motifs_order <- order.dendrogram(clustering_motifs)
  
  dendro_motifs <- clustering_motifs %>%
    dendro_data %>%
    segment %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_segment() +
    scale_y_continuous("") +
    scale_x_discrete("",
                     limits = factor(1L:nobs(clustering_motifs))) + 
    theme_void() + 
    coord_flip()
  
  test_res[["motif"]] <- factor(test_res[["motif"]],
                                levels = test_res_plot_dat[["motif"]][motifs_order], ordered = TRUE)
  
  clustering_comparisons <- as.dendrogram(hclust(dist(t(as.matrix(test_res_plot_dat[,2:7])))))
  comparisons_order <- order.dendrogram(clustering_comparisons)
  dendro_comparisons <- clustering_comparisons %>%
    dendro_data %>%
    segment %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_segment() +
    scale_y_continuous("") +
    scale_x_discrete("",
                     limits = factor(1L:nobs(clustering_comparisons))) + 
    theme_void() 
  
  test_res[["comparison"]] <- factor(test_res[["comparison"]],
                                     levels = colnames(test_res_plot_dat)[2:7][comparisons_order], ordered = TRUE)
  
  heatmap <- test_res %>% 
    select(-pval)  %>% 
    group_by(motif, comparison) %>% 
    summarise(is_significant = as.logical(pval_adjusted < 0.05)) %>% 
    ggplot(aes(y = motif, x = comparison)) +
    geom_tile(aes(fill = is_significant)) +
    scale_fill_manual(values = c("TRUE" = "red4", "FALSE" = "wheat")) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.key.width = unit(2, "lines"),
          axis.text.y = element_text(size = 3)) 
  
  max_height <- unit.pmax(ggplotGrob(heatmap)[["heights"]],
                          ggplotGrob(dendro_motifs)[["heights"]])
  
  grob_list <- list(heatmap = ggplotGrob(heatmap),
                    dendrogram_right = ggplotGrob(dendro_motifs),
                    dendrogram_top = ggplotGrob(dendro_comparisons))
  
  max_width <- unit.pmax(grob_list[["heatmap"]][["widths"]],
                         grob_list[["dendrogram_top"]][["widths"]])
  
  grob_list[["heatmap"]][["widths"]] <-
    grob_list[["dendrogram_top"]][["widths"]] <-
    max_width
  
  grob_list[["heatmap"]][["heights"]] <-
    grob_list[["dendrogram_right"]][["heights"]] <-
    max_height
  
  top_right <- grid.rect(gp = gpar(col = NA), draw = FALSE)
  
  grid.arrange(grob_list[["dendrogram_top"]],
               top_right,
               grob_list[["heatmap"]],
               grob_list[["dendrogram_right"]],
               widths = c(0.9, 0.1), heights = c(0.025, 0.975))
}
