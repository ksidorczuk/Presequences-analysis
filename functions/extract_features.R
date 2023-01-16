#' Calculate ngram presence/absence matrix
#' 
#' This function calculates a presence/absence matrix
#' of ngrams with given length and gaps for provided
#' sequences.
#' @param sequences a list of sequences to analyse
#' @param k_vector a vector of ngram lengths to calculate
#' @param kmer_gaps_list a list of gaps in ngrams to use
#' @return a data frame of ngram binary presence/absence 
#' with ngrams in columns and sequences in rows
calculate_binary_ngram_matrix <- function(sequences, k_vector, kmer_gaps_list) {
  count_multimers(sequences,
                  k_vector = k_vector, # = c(1, 2, 2, rep(3, 4)),
                  kmer_gaps_list = kmer_gaps_list, # = list(NULL, NULL, 1, c(0, 0), c(1, 0), c(0, 1), c(1, 1)),
                  kmer_alphabet = toupper(colnames(biogram::aaprop)),
                  with_kmer_counts = FALSE) %>%
    as.matrix() %>% 
    as.data.frame()
}

#' Encode a sequence
#' 
#' Encode a sequence using a physicochemical property index code
#' from AAindex database.
#' @param x a list of sequences to calculate the values
#' @param property a name of AAindex physicochemical property
#' to use for encoding, e.g. \code{"BIGC670101"} for Residue volume
#' (Bigelow, 1967).
#' @return mean values of given property for each analysed sequence
encode_seq <- function(x, property) {
  sapply(x, function(ith_seq) {
    mean(aaprop[property, tolower(ith_seq)])
  })
}

#' Calculate physicochemical properties
#' 
#' This function calculates physicochemical properties
#' for many properties using \code{\link{encode_seq}}
#' function.
#' @param sequences sequences for which the properties
#' are to be calculated
#' @param prop_list a list of property names. AAindex codes
#' should be used, e.g. \code{"BIGC670101"} for Residue volume
#' (Bigelow, 1967).
#' @return a data frame of physicochemical property
#' values calculated for sequences, where each row corresponds
#' to one sequence and columns correspond to properties.
calculate_properties <- function(sequences, prop_list) {
  lapply(prop_list, function(ith_prop) {
    encode_seq(sequences, ith_prop) 
  }) %>% setNames(prop_list) %>% do.call(cbind, .)
}

#' Calculate amino acid composition of peptides
#' 
#' This function calculates amino acid composition
#' of all peptides from a list of datasets.
#' @param dataset_list a list of datasets for which
#' amino acid composition will be calculated
#' @return a data frame with two columns corresponding
#' amino acid and its frequency, where each row corresponds
#' to a peptide.
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
