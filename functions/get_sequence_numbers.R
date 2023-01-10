#' Summarise transit peptide information
#' 
#' This function summarizes information about sequences with transit peptides.
#' Based on the names of lists of sequences and their annotations from UniProt 
#' it analyzes these lists and returns information about number of annotated
#' transit peptides, transit peptides containing nonstandard amino acids, 
#' shorter than 10 and the number of sequences to use in further analyses
#' (only standard amino acids, length >= 10).
#' 
#' @param seq_dataset_vec vector of names of sequence datasets in a form
#' of a list 
#' @param annotations a data frame of annotations from UniProt
#' @return a data frame with five columns describing: dataset name,
#' number of sequences containing transit peptides, number of transit
#' peptides with nonstandard amino acids, number of transit peptides
#' shorter than 10, number of transit peptides suitable for further
#' analyses. Each row corresponds to one dataset.
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

#' Summarise signal peptide information
#' 
#' This function summarizes information about sequences with signal peptides.
#' Based on the names of lists of sequences and their annotations from UniProt 
#' it analyzes these lists and returns information about number of annotated
#' signal peptides, signal peptides containing nonstandard amino acids, 
#' shorter than 10 and the number of sequences to use in further analyses
#' (only standard amino acids, length >= 10).
#' 
#' @param seq_dataset_vec vector of names of sequence datasets in a form
#' of a list 
#' @param annotations a data frame of annotations from UniProt
#' @return a data frame with five columns describing: dataset name,
#' number of sequences containing signal peptides, number of signal
#' peptides with nonstandard amino acids, number of signal peptides
#' shorter than 10, number of signal peptides suitable for further
#' analyses. Each row corresponds to one dataset.
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

#' Summarise transmembrane domain information
#' 
#' This function summarizes information about sequences with transmembrane domains.
#' Based on the names of lists of sequences and their annotations from UniProt 
#' it analyzes these lists and returns information about number of annotated
#' transmembrane regions with alpha helical or beta sheet structure, how many
#' of them contain nonstandard amino acids, are shorter than 10 and how many
#' are suitable for further analyses (only standard amino acids, length >= 10).
#' 
#' @param seq_dataset_vec vector of names of sequence datasets in a form
#' of a list 
#' @param annotations a data frame of annotations from UniProt
#' @return a data frame with five columns describing: dataset name,
#' number of sequences containing helical transmembrane domains,
#' containing beta transmembrane domains, helical containing nonstandard
#' amino acids, beta containing nonstandard amino acids, helical shorter than 10,
#' beta shorter than 10, helical suitable for analyses and beta suitable for
#' analyses. Each row corresponds to one dataset.
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

#' Summarise amyloid domain information
#' 
#' This function summarizes information about amyloid datasets.
#' Based on the named list of sequence datasets it analyzes them and returns 
#' information about number of sequences, sequences with nonstandard amino acids,
#' shorter than 10 and suitable for further analyses (only standard amino acids, 
#' length >= 10).
#' 
#' @param dataset_list named list of sequence datasets, where each dataset
#' is a sequence list
#' @return a data frame with five columns describing: dataset name,
#' number of sequences, containing nonstandard amino acids, shorter than 10
#' and suitable for analyses. Each row corresponds to one dataset.
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

#' Summarise AMP information
#' 
#' This function summarizes information about AMP dataset.
#' Based on the sequence datasets it analyzes it and returns information 
#' about number of sequences, sequences with nonstandard amino acids,
#' shorter than 10 and suitable for further analyses (only standard amino acids, 
#' length >= 10).
#' 
#' @param amp_all sequence dataset in a form of a list
#' @return a data frame with six columns describing: dataset name,
#' number of sequences, containing nonstandard amino acids, shorter 
#' than 10, longer than 10 and suitable for analyses.
summarise_amp_info <- function(amp_all) {
  data.frame(
    dataset = "AMP",
    sequences = length(amp_all),
    with_nonstandard_aa = length(amp_all) - length(filter_nonstandard_aa(amp_all)),
    shorter_than_10 = sum(lengths(amp_all) < 10),
    longer_than_100 = sum(lengths(amp_all) > 100),
    to_analyse = sum(lengths(filter_nonstandard_aa(amp_all)) >= 10)
  )
}
