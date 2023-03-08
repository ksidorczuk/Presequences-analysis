#' Count the longest stretch of AMP mers
#' 
#' Calculates the length of the longest stretch of
#' consecutive mers predicted as possessing antimicrobial
#' properties
#' @param x a numeric vector of prediction values
#' @return length of the longest stretch of mers predicted as AMPs
count_longest <- function(x) {
  splitted_x <- strsplit(x = paste0(as.numeric(x > 0.5), collapse = ""),
                         split = "0")[[1]]
  len <- unname(sapply(splitted_x, nchar))
  if (length(len[len > 0]) == 0) {
    0 } else {
      len[len > 0]
    }
}

#' Calculate AmpGram statistics
#' 
#' Calculates statistics from mer predictions required 
#' for the second layer of the AmpGram model. 
#' @param mer_preds AmpGram predictions for mets
#' @return a data frame of statistics calculated for
#' each analysed peptide
calculate_statistics <- function(mer_preds) {
  group_by(mer_preds, source_peptide) %>% 
    summarise(fraction_true = mean(pred > 0.5),
              pred_mean = mean(pred),
              pred_median = median(pred),
              n_peptide = length(pred),
              n_pos = sum(pred > 0.5),
              pred_min = min(pred),
              pred_max = max(pred), 
              longest_pos = max(count_longest(pred)),
              n_pos_10 = sum(count_longest(pred) >= 10),
              frac_0_0.2 = sum(pred <= 0.2)/length(pred),
              frac_0.2_0.4 = sum(pred > 0.2 & pred <= 0.4)/length(pred),
              frac_0.4_0.6 = sum(pred > 0.4 & pred <= 0.6)/length(pred),
              frac_0.6_0.8 = sum(pred > 0.6 & pred <= 0.8)/length(pred),
              frac_0.8_1 = sum(pred > 0.8 & pred <= 1)/length(pred)) 
}

#' Predict AMPs using AmpGram
#' 
#' This function predicts antimicrobial properties of peptides
#' listed in a FASTA-formatted file and saves the results to specified
#' output directory and file. It runs on chunks of the sequence file
#' analysing up to 20,000 sequences at once to allow analysis of very 
#' large datasets step by step. All results are then combined. 
#' @param sequence_file name of a file containing sequences in a FASTA
#' format to analyse
#' @param output_path path of the output directory
#' @param output_name name of the output file
predict_with_AmpGram <- function(sequence_file, output_path, output_name) {
  seqs <- read_fasta(sequence_file)
  lapply(seq(1, length(seqs), 20000), function(i) {
    if(i+19999 < length(seqs)) {
      x <- i+19999
    } else {
      x <- length(seqs)
    }
    dat <- seqs[i:x]
    mers <- lapply(seq_along(dat), function(ith_seq) {
      lapply(1:(length(dat[[ith_seq]])-9), function(i) dat[[ith_seq]][i:(i+9)]) %>% 
        setNames(rep(names(dat)[ith_seq], length(dat[[ith_seq]])-9))
    }) %>% unlist(recursive = FALSE)
    ngrams <- count_multimers(mers,
                              k_vector = c(1, rep(2, 4), c(rep(3, 4))),
                              kmer_gaps_list = list(NULL, NULL, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
                              kmer_alphabet = toupper(colnames(aaprop)),
                              with_kmer_counts = FALSE)
    colnames(ngrams)[1:20] <- paste0(colnames(ngrams)[1:20], "_0")
    to_add <- AmpGram_model[["imp_features"]][which(!(AmpGram_model[["imp_features"]] %in% colnames(ngrams)))]
    ngrams <- cbind(ngrams, matrix(0, ncol = length(to_add), nrow = nrow(ngrams), dimnames = list(NULL, to_add)))
    all_mers_pred <- data.frame(source_peptide = names(mers),
                                pred = predict(AmpGram_model[["rf_mers"]], 
                                               ngrams[, AmpGram_model[["imp_features"]]])[["predictions"]][, 2])
    single_prot_pred <- data.frame(source_peptide = unique(all_mers_pred[["source_peptide"]]),
                                   pred = predict(AmpGram_model[["rf_peptides"]], 
                                                  calculate_statistics(all_mers_pred))[["predictions"]][, 2]) %>% 
      mutate(decision = ifelse(pred > 0.5, TRUE, FALSE))
    write.csv(single_prot_pred, file = paste0(output_path, "AmpGram_", i, "_temp.csv"), row.names = FALSE)
  })
  res <- list.files(output_path, pattern = "_temp.csv", full.names = TRUE)
  all_res <- lapply(res, function(ith_file) {
    read.csv(ith_file)
  }) %>% bind_rows()
  file.remove(res)
  write.csv(all_res, output_name, row.names = FALSE)
}
