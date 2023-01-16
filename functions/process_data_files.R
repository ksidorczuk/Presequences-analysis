#' Read AmyPro data
#' 
#' This function reads in the raw dataset from AmyPro database and transforms
#' it into a data frame.
#' 
#' @param amypro_txt_file path to the raw txt file downloaded from the AmyPro
#' database
#' @return a data frame with AmyPro sequences and annotations
read_amypro_dat <- function(amypro_txt_file) {
  read.table(amypro_txt_file, sep = "\t", fill = TRUE, row.names = NULL) %>% 
    setNames(colnames(.)[2:14]) %>% 
    .[,1:13]
}

#' Extract AmyPro regions
#' 
#' This function extracts amyloidogenic regions from proteins from AmyPro database.
#' 
#' @param amypro_dat a data frame with AmyPro data
#' @return a list of amyloidogenic regions with names corresponding to the 
#' original protein from which it was extracted.
get_amypro_regions <- function(amypro_dat) {
  lapply(1:nrow(amypro_dat), function(i) {
    if(amypro_dat[["regions"]][i] != "-") {
      x <- strsplit(amypro_dat[["regions"]][i], ",")[[1]]
      region_list <- as.list(setNames(x, paste0(amypro_dat[['ID']][i], "_", 1:length(x))))
      lapply(region_list, function(ith_region) strsplit(ith_region, "")[[1]])
    }
  }) %>% unlist(recursive = FALSE)
}

#' Extract CPAD peptides
#' 
#' This function extracts peptides classified as amyloids in the CPAD database.
#' @param cpad_dat a data frame of data downloaded from the CPAD database
#' @return a list of amyloid peptides with names corresponding to CPAD entry names.
get_cpad_peptides <- function(cpad_dat) {
  amyloids <- cpad_dat %>% 
    filter(Classification == "Amyloid")
  amyloids[["Peptide"]] %>% 
    setNames(amyloids[["Entry"]]) %>% 
    lapply(., function(ith_pep) strsplit(ith_pep, "")[[1]])
}

#' Extract AMPs
#' 
#' This function extracts sequences of AMPs from the DBAASP dataset. 
#' @param dbaasp_dat a data frame of data downloaded from the DBAASP database
#' @return a list of AMP sequences with the names corresponding to
#' DBAASP ids. 
get_AMPs <- function(dbaasp_dat) {
  dbaasp_dat[["sequence"]] %>% 
    lapply(., function(ith_pep) strsplit(ith_pep, "")[[1]]) %>% 
    setNames(dbaasp_dat[["dbaasp_id"]])
}

#' Filter presequence entries
#' 
#' Filter UniProt annotations to select only entries corresponding
#' to presequences in the datasets
#' @param datasets_list list of named datasets
#' @param annotations_all a data frame of all UniProt annotations
#' @return a data frame of UniProt annotations concerning presequences
filter_presequence_entries <- function(datasets_list, annotations_all) {
  presequences_acc <- datasets_list[c("cTP experimentally verified presequence", "mTP experimentally verified presequence",
                                      "cTP-mTP experimentally verified presequence", "SP experimentally verified presequence")] %>% 
    unname() %>% 
    unlist(recursive = FALSE) %>% 
    names()
  
  annotations_all %>% 
    filter(Entry %in% presequences_acc) 
}

