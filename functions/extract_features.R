calculate_binary_ngram_matrix <- function(sequences, k_vector, kmer_gaps_list) {
  count_multimers(sequences,
                  k_vector = k_vector, # = c(1, 2, 2, rep(3, 4)),
                  kmer_gaps_list = kmer_gaps_list, # = list(NULL, NULL, 1, c(0, 0), c(1, 0), c(0, 1), c(1, 1)),
                  kmer_alphabet = toupper(colnames(biogram::aaprop)),
                  with_kmer_counts = FALSE) %>%
    as.matrix() %>% 
    as.data.frame()
}

encode_seq <- function(x, property) {
  sapply(x, function(ith_seq) {
    mean(aaprop[property, tolower(ith_seq)])
  })
}

calculate_properties <- function(sequences, prop_list) {
  lapply(prop_list, function(ith_prop) {
    encode_seq(sequences, ith_prop) 
  }) %>% setNames(prop_list) %>% do.call(cbind, .)
}

calculate_aa_comp_peptides <- function(dataset_list) {
  lapply(names(dataset_list), function(i) {
    lapply(names(dataset_list[[i]]), function(ith_prot) {
      data.frame(table(factor(dataset_list[[i]][[ith_prot]], levels = c("A", "C", "D", "E", "F", "G", "H",
                                                                        "I", "K", "L", "M", "N", "P", "Q",
                                                                        "R", "S", "T", "V", "W", "Y")))/length(dataset_list[[i]][[ith_prot]])) %>%
        setNames(c("Amino acid", "Frequency")) %>% 
        mutate(dataset = i,
               prot = ith_prot)
    }) %>% bind_rows()
  }) %>% bind_rows()
} 
