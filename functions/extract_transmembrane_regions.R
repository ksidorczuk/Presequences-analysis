extract_transmembrane_regions <- function(sequences, annotation_df, remove_nonstandard = TRUE) {
  
  names(sequences) <- sapply(names(sequences), function(i) strsplit(i, "|", fixed = TRUE)[[1]][2])
  
  dat <- annotation_df %>% 
    filter(!(is.na(`Transmembrane`))) %>% 
    select(c(Entry, `Transmembrane`)) %>% 
    filter(Entry %in% names(sequences))
  
  tm_dat <- lapply(dat[["Entry"]], function(ith_acc) {
    x <- strsplit(dat[["Transmembrane"]][which(dat[["Entry"]] == ith_acc)], "TRANSMEM ", fixed = TRUE)[[1]][-1] %>% 
      .[which(grepl("Helical|Beta stranded", .) & (!(grepl("<1|>|-2:", .))) & (!grepl("?", ., fixed = TRUE)))]

    if(length(x) > 0) {
      lapply(1:length(x), function(i) {
        data.frame(seq = ith_acc,
                   tm_name = paste0(ith_acc, "_TM_", i),
                   tm_start = strsplit(strsplit(x[i], ";")[[1]][1], "..", fixed = TRUE)[[1]][1],
                   tm_end = strsplit(strsplit(x[i], ";")[[1]][1], "..", fixed = TRUE)[[1]][2],
                   type = ifelse(grepl("Helical", x[i]), "Helical", "Beta stranded"))
      }) %>% bind_rows()
    }
  }) %>% bind_rows()
  tm_dat_a <- filter(tm_dat, type == "Helical")
  tm_dat_b <- filter(tm_dat, type == "Beta stranded")
  
  tms_a <- lapply(1:nrow(tm_dat_a), function(j) {
    sequences[[tm_dat_a[["seq"]][j]]][(tm_dat_a[["tm_start"]][j]):(tm_dat_a[["tm_end"]][j])]
  }) %>% setNames(tm_dat_a[["tm_name"]])
  
  tms_b <- lapply(1:nrow(tm_dat_b), function(j) {
    sequences[[tm_dat_b[["seq"]][j]]][(tm_dat_b[["tm_start"]][j]):(tm_dat_b[["tm_end"]][j])]
  }) %>% setNames(tm_dat_b[["tm_name"]])
  
  if(isTRUE(remove_nonstandard)) {
    tms_a <- filter_nonstandard_aa(tms_a)
    tms_b <- filter_nonstandard_aa(tms_b)
  }
  
  list("alpha" = tms_a,
       "beta" = tms_b)
}


extract_intramembrane_regions <- function(sequences, annotation_df, remove_nonstandard = TRUE) {
  
  names(sequences) <- sapply(names(sequences), function(i) strsplit(i, "|", fixed = TRUE)[[1]][2])
  
  dat <- annotation_df %>% 
    filter(!(is.na(`Intramembrane`))) %>% 
    select(c(Entry, `Intramembrane`)) %>% 
    filter(Entry %in% names(sequences))
  
  intm_dat <- lapply(dat[["Entry"]], function(ith_acc) {
    x <- strsplit(dat[["Intramembrane"]][which(dat[["Entry"]] == ith_acc)], "INTRAMEM ", fixed = TRUE)[[1]][-1] %>% 
      .[which(grepl("Helical|Pore-forming", .) & (!(grepl("?..", ., fixed = TRUE))))] 
    
    lapply(1:length(x), function(i) {
      data.frame(seq = ith_acc,
                 tm_name = paste0(ith_acc, "_INTM_", i),
                 tm_start = strsplit(strsplit(x[i], ";")[[1]][1], "..", fixed = TRUE)[[1]][1],
                 tm_end = strsplit(strsplit(x[i], ";")[[1]][1], "..", fixed = TRUE)[[1]][2],
                 type = ifelse(grepl("Helical", x[i]), "Helical", "Pore-forming"))
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  intm_dat_a <- filter(intm_dat, type == "Helical")
  intm_dat_p <- filter(intm_dat, type == "Pore-forming")
  
  intms_a <- lapply(1:nrow(intm_dat_a), function(j) {
    sequences[[intm_dat_a[["seq"]][j]]][(intm_dat_a[["tm_start"]][j]):(intm_dat_a[["tm_end"]][j])]
  }) %>% setNames(tm_dat_a[["tm_name"]])
  
  intms_p <- lapply(1:nrow(tm_dat_b), function(j) {
    sequences[[intm_dat_p[["seq"]][j]]][(intm_dat_p[["tm_start"]][j]):(intm_dat_p[["tm_end"]][j])]
  }) %>% setNames(tm_dat_b[["tm_name"]])
  
  if(isTRUE(remove_nonstandard)) {
    intms_a <- filter_nonstandard_aa(intms_a)
    intms_p <- filter_nonstandard_aa(intms_p)
  }
  
  list("helical" = intms_a,
       "pore-forming" = intms_p)
}
