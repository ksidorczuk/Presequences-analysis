calculate_binary_ngram_matrix <- function(sequences) {
  count_multimers(sequences,
                  k_vector = c(1, 2, 2, rep(3, 4)),
                  kmer_gaps_list = list(NULL, NULL, 1, c(0, 0), c(1, 0), c(0, 1), c(1, 1)),
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
