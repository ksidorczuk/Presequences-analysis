#' Filter sequences with nonstandard amino acids
#' 
#' Analyzes a dataset of sequences and filters out those
#' containing amino acids other than standard.
#' 
#' @param sequences a list of amino acid sequences
#' @return a list of amino acid sequences containing only
#' standard amino acids
filter_nonstandard_aa <- function(sequences) {
  standard <- toupper(biogram:::return_elements(seq_type = "prot"))
  is_standard <- vapply(sequences, function(seq) all(seq %in% standard), c(FALSE))
  sequences[is_standard]
}


#' Extract transit peptides
#' 
#' This function extracts transit peptides from sequences based on UniProt
#' annotations. It considers only the transit peptides with properly annotated
#' start and end positions and not truncated.
#' 
#' @param sequences a list of sequences
#' @param annotation_df a data frame of UniProt annotations
#' @param remove_nonstandard a logical indicating if transit peptides containing
#' nonstandard amino acids should be removed
#' @return a list of extracted transit peptides with names corresponding to their
#' original sequences.
extract_transit_peptides <- function(sequences, annotation_df, remove_nonstandard = TRUE) {
  
 names(sequences) <- sapply(names(sequences), function(i) strsplit(i, "|", fixed = TRUE)[[1]][2])
  
  dat <- annotation_df %>% 
    filter(!is.na(`Transit peptide`)) %>% 
    select(c(Entry, `Transit peptide`)) %>% 
    filter(Entry %in% names(sequences)  & (!(grepl("1..?", `Transit peptide`, fixed = TRUE) | grepl("2..?", `Transit peptide`, fixed = TRUE))) & (!(grepl("<1", `Transit peptide`))) & (!(grepl("TRANSIT 1;", `Transit peptide`))))
  
  position_list <- lapply(dat[["Entry"]], function(ith_acc) {
    x <- strsplit(dat[["Transit peptide"]][which(dat[["Entry"]] == ith_acc)], ".", fixed = TRUE)[[1]][3]
    x <- strsplit(x, ";")[[1]][1]
    setNames(x, ith_acc)
  }) %>% unlist()
  
  tps <- lapply(names(position_list), function(ith_prot) {
    sequences[[ith_prot]][1:position_list[[ith_prot]]]
  }) %>% setNames(names(position_list))
  
  if(isTRUE(remove_nonstandard)) {
    tps <- filter_nonstandard_aa(tps)
  }
  
  tps
}


#' Extract signal peptides
#' 
#' This function extracts signal peptides from sequences based on UniProt
#' annotations. It considers only the signal peptides with properly annotated
#' start and end positions and not truncated.
#' 
#' @param sequences a list of sequences
#' @param annotation_df a data frame of UniProt annotations
#' @param remove_nonstandard a logical indicating if signal peptides containing
#' nonstandard amino acids should be removed
#' @return a list of extracted signal peptides with names corresponding to their
#' original sequences.
extract_signal_peptides <- function(sequences, annotation_df, remove_nonstandard = TRUE) {
  
  names(sequences) <- sapply(names(sequences), function(i) strsplit(i, "|", fixed = TRUE)[[1]][2])
  
  dat <- annotation_df %>% 
    filter(!(is.na(`Signal peptide`))) %>% 
    select(c(Entry, `Signal peptide`)) %>% 
    filter(Entry %in% names(sequences)  & (!(grepl("1..?", `Signal peptide`, fixed = TRUE) | grepl("2..?", `Signal peptide`, fixed = TRUE) | grepl("<1", `Signal peptide`) | grepl("SIGNAL 1;", `Signal peptide`) | grepl("..>", `Signal peptide`, fixed = TRUE))))
  
  position_list <- lapply(dat[["Entry"]], function(ith_acc) {
    x <- strsplit(dat[["Signal peptide"]][which(dat[["Entry"]] == ith_acc)], ".", fixed = TRUE)[[1]][3]
    x <- strsplit(x, ";")[[1]][1]
    setNames(x, ith_acc)
  }) %>% unlist()
  
  sps <- lapply(names(position_list), function(ith_prot) {
    sequences[[ith_prot]][1:position_list[[ith_prot]]]
  }) %>% setNames(names(position_list))
  
  if(isTRUE(remove_nonstandard)) {
    sps <- filter_nonstandard_aa(sps)
  }
  
  sps
}

