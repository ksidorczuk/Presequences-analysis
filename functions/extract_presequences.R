filter_nonstandard_aa <- function(sequences) {
  standard <- toupper(biogram:::return_elements(seq_type = "prot"))
  is_standard <- vapply(sequences, function(seq) all(seq %in% standard), c(FALSE))
  sequences[is_standard]
}


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

