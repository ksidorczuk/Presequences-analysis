read_amypro_dat <- function(amypro_txt_file) {
  read.table(amypro_txt_file, sep = "\t", fill = TRUE, row.names = NULL) %>% 
    setNames(colnames(.)[2:14]) %>% 
    .[,1:13]
}

get_amypro_regions <- function(amypro_dat) {
  lapply(1:nrow(amypro_dat), function(i) {
    if(amypro_dat[["regions"]][i] != "-") {
      x <- strsplit(amypro_dat[["regions"]][i], ",")[[1]]
      region_list <- as.list(setNames(x, paste0(amypro_dat[['ID']][i], "_", 1:length(x))))
      lapply(region_list, function(ith_region) strsplit(ith_region, "")[[1]])
    }
  }) %>% unlist(recursive = FALSE)
}

get_cpad_peptides <- function(cpad_dat) {
  amyloids <- cpad_dat %>% 
    filter(Classification == "Amyloid")
  amyloids[["Peptide"]] %>% 
    setNames(amyloids[["Entry"]]) %>% 
    lapply(., function(ith_pep) strsplit(ith_pep, "")[[1]])
}

get_AMPs <- function(dbaasp_dat) {
  dbaasp_dat[["sequence"]] %>% 
    lapply(., function(ith_pep) strsplit(ith_pep, "")[[1]]) %>% 
    setNames(dbaasp_dat[["dbaasp_id"]])
}


filter_presequence_entries <- function(datasets_list, annotations_all) {
  presequences_acc <- datasets_list[c("cTP experimentally verified presequence", "mTP experimentally verified presequence",
                                      "cTP-mTP experimentally verified presequence", "SP experimentally verified presequence")] %>% 
    unname() %>% 
    unlist(recursive = FALSE) %>% 
    names()
  
  annotations_all %>% 
    filter(Entry %in% presequences_acc) 
}

