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
calculate_ngram_freq <- function(datasets, k, kmer_gaps) {
  freqs <- lapply(names(datasets), function(ith_dataset) {
    x <- calculate_binary_ngram_matrix(datasets[[ith_dataset]],
                                       k = k,
                                       kmer_gaps_list = kmer_gaps) 
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
#' @param datasets a named list of sequence datasets 
#' @param motifs a vector of motifs to find (written as
#' regular expressions)
#' 
#' @return a data frame with the following columns: start/end 
#' positions of a motif, dataset name, sequence name and motif
locate_motifs <- function(datasets, motifs) {
  pblapply(motifs, cl = 8, function(ith_motif) {
    lapply(names(datasets), function(ith_set) {
      lapply(names(datasets[[ith_set]]), function(ith_seq) {
        str_locate_all(paste0(datasets[[ith_set]][[ith_seq]], collapse = ""), ith_motif) %>%
          as.data.frame() %>%
          mutate(dataset = ith_set,
                 seq = ith_seq,
                 motif = ith_motif)
      }) %>% bind_rows()
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
get_motif_plot_dat <- function(motif_freqs, gapped = FALSE) {
  if(gapped == TRUE) {
    colnames(motif_freqs) <- sapply(colnames(motif_freqs), function(i) {
      ifelse(i == "dataset", "dataset", gsub("_", ".", decode_ngrams(i), fixed = TRUE))
    })
  } else {
    colnames(motif_freqs) <- sapply(colnames(motif_freqs), function(i) {
      ifelse(i == "dataset", "dataset", gsub(".", "", gsub("_0", "", i), fixed = TRUE))
    }) 
  }
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
  
  plot_dat <- pblapply(1:nrow(lens_dat), cl = 8, function(i) {
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
  
  lens_dat <- left_join(filter(found_motifs, dataset == dataset_name), lens, by = "seq")
  
  plot_dat <- pblapply(1:nrow(lens_dat), cl = 8, function(i) {
    data.frame(dataset = lens_dat[["dataset"]][i],
               seq_name = lens_dat[["seq"]][i],
               motif_pos = lens_dat[["start"]][i]:lens_dat[["end"]][i],
               len = lens_dat[["len"]][i],
               motif = lens_dat[["motif"]][[i]])
  }) %>% bind_rows()
  
  n_ngrams_plot_dat <- plot_dat %>% 
    group_by(dataset, seq_name, len, motif, motif_pos) %>% 
    summarise(n = n()) 
  
  pblapply(1:nrow(lens_dat), cl = 8, function(i) {
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


plot_motif_positions_scaled <- function(dataset_list, dataset_name, found_motifs, type = "point") {
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
      geom_violin() +
      facet_wrap(~motif) +
      theme_bw()  +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      ggtitle(dataset_name)
  }

}
