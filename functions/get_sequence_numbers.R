summarise_TP_info <- function(seq_dataset_vec, annotations) {
  lapply(seq_dataset_vec, function(i) {
    tps <- extract_transit_peptides(sequences = get(i),
                                    annotation_df = annotations,
                                    remove_nonstandard = FALSE)
    standard <- filter_nonstandard_aa(tps)
    data.frame(
      dataset = i, 
      TP = length(tps),
      with_nonstandard_aa = length(tps) - length(standard),
      shorter_than_10 = sum(lengths(tps) < 10),
      to_analyse = sum(lengths(standard) >= 10)
    )
  }) %>% bind_rows()
}

summarise_SP_info <- function(seq_dataset_vec, annotations) {
  lapply(seq_dataset_vec, function(i) {
    sps <- extract_signal_peptides(sequences = get(i),
                                   annotation_df = annotations,
                                   remove_nonstandard = FALSE)
    standard <- filter_nonstandard_aa(sps)
    data.frame(
      dataset = i, 
      SP = length(sps),
      with_nonstandard_aa = length(sps) - length(standard),
      shorter_than_10 = sum(lengths(sps) < 10),
      to_analyse = sum(lengths(standard) >= 10)
    )
  }) %>% bind_rows()
}

summarise_TM_info <- function(seq_dataset_vec, annotations) {
  lapply(seq_dataset_vec, function(i) {
    tms <- extract_transmembrane_regions(sequences = get(i),
                                         annotation_df = annotations,
                                         remove_nonstandard = FALSE)
    standard_a <- filter_nonstandard_aa(tms[["alpha"]])
    standard_b <- filter_nonstandard_aa(tms[["beta"]])
    data.frame(
      dataset = i, 
      TM_helical = length(tms[["alpha"]]),
      TM_beta = length(tms[["beta"]]),
      with_nonstandard_aa_helical = length(tms[["alpha"]]) - length(standard_a),
      with_nonstandard_aa_beta = length(tms[["beta"]]) - length(standard_b),
      shorter_than_10_helical = sum(lengths(tms[["alpha"]]) < 10),
      shorter_than_10_beta = sum(lengths(tms[["beta"]]) < 10),
      to_analyse_helical = sum(lengths(standard_a) >= 10),
      to_analyse_beta = sum(lengths(standard_b) >= 10)
    )
  }) %>% bind_rows()
}

summarise_amyloid_info <- function(dataset_list) {
  lapply(names(dataset_list), function(ith_set) {
    data.frame(
      dataset = ith_set,
      sequences = length(dataset_list[[ith_set]]),
      with_nonstandard_aa = length(dataset_list[[ith_set]]) - length(filter_nonstandard_aa(dataset_list[[ith_set]])),
      shorter_than_10 = sum(lengths(dataset_list[[ith_set]]) < 10),
      to_analyse = sum(lengths(dataset_list[[ith_set]]) >= 10)
    )
  }) %>% bind_rows()
}

summarise_amp_info <- function(amp_all) {
  data.frame(
    dataset = "AMP",
    sequences = length(amp_all),
    with_nonstandard_aa = length(amp_all) - length(filter_nonstandard_aa(amp_all)),
    shorter_than_10 = sum(lengths(amp_all) < 10),
    to_analyse = sum(lengths(filter_nonstandard_aa(amp_all)) >= 10)
  )
}
