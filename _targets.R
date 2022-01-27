library(targets)
library(dplyr)
library(readxl)
library(biogram)

data_path <- "/media/kasia/Data/Dropbox/Presequences/"

source("functions/extract_presequences.R")
source("functions/extract_transmembrane_regions.R")
source("functions/get_sequence_numbers.R")
source("functions/process_data_files.R")
source("functions/filter_sequences.R")

list(
  # Annotation files
  tar_target(
    annotations_all_file,
    paste0(data_path, "Raw_data/uniprot-reviewed yes.xlsx"),
    format = "file"
  ),
  tar_target(
    annotations_tm_file,
    paste0(data_path, "Raw_data/uniprot-annotation (type transmem)-filtered-reviewed yes.xlsx"),
    format = "file"
  ),
  tar_target(
    annotations_intm_file,
    paste0(data_path, "Raw_data/uniprot-annotation (type intramem)+reviewed yes.xlsx"),
    format = "file"
  ),
  # Annotations
  tar_target(
    annotations_all,
    read_xlsx(annotations_all_file)
  ),
  tar_target(
    annotations_tm,
    read_xlsx(annotations_tm_file)
  ),
  # Sequence files
  tar_target(
    `cTP_loc_tp.exp_file`,
    paste0(data_path, "Raw_data/cTP_loc_tp-exp.fasta"),
    format = "file"
  ),
  tar_target(
    cTP_loc.exp_tp_file,
    paste0(data_path, "Raw_data/cTP_loc-exp_tp.fasta"),
    format = "file"
  ),
  tar_target(
    `cTP.mTP_loc_tp.exp_file`,
    paste0(data_path, "Raw_data/cTP-mTP_loc_tp-exp.fasta"),
    format = "file"
  ),
  tar_target(
    `cTP.mTP_loc.exp_tp_file`,
    paste0(data_path, "Raw_data/cTP-mTP_loc-exp_tp.fasta"),
    format = "file"
  ),
  tar_target(
    `mTP_loc.exp_tp_file`,
    paste0(data_path, "Raw_data/mTP_loc-exp_tp.fasta"),
    format = "file"
  ),
  tar_target(
    `mTP_loc_tp.exp_file`,
    paste0(data_path, "Raw_data/mTP_loc_tp-exp.fasta"),
    format = "file"
  ),
  tar_target(
    `SP_secreted_sp.exp_file`,
    paste0(data_path, "Raw_data/SP_secreted_sp-exp.fasta"),
    format = "file"
  ),
  tar_target(
    `SP_secreted.er_sp.exp_file`,
    paste0(data_path, "Raw_data/SP_secreted-er_sp-exp.fasta"),
    format = "file"
  ),
  tar_target(
    `SP_secreted.er.exp_sp_file`,
    paste0(data_path, "Raw_data/SP_secreted-er-exp_sp.fasta"),
    format = "file"
  ),
  tar_target(
    `SP_secreted.exp_sp_file`,
    paste0(data_path, "Raw_data/SP_secreted-exp_sp.fasta"),
    format = "file"
  ),
  tar_target(
    `SP_sp.exp_file`,
    paste0(data_path, "Raw_data/SP_sp-exp.fasta"),
    format = "file"
  ),
  tar_target(
    `SP_sp_file`,
    paste0(data_path, "Raw_data/SP_sp.fasta"),
    format = "file"
  ),
  tar_target(
    `TM_exp_file`,
    paste0(data_path, "Raw_data/TM_exp.fasta"),
    format = "file"
  ),
  tar_target(
    `TM_file`,
    paste0(data_path, "Raw_data/TM.fasta"),
    format = "file"
  ),
  # Sequences files
  tar_target(
    cTP_loc_tp.exp_seqs,
    read_fasta(cTP_loc_tp.exp_file)
  ),
  tar_target(
    cTP_loc.exp_tp_seqs,
    read_fasta(cTP_loc.exp_tp_file)
  ),
  tar_target(
    cTP.mTP_loc_tp.exp_seqs,
    read_fasta(cTP.mTP_loc_tp.exp_file)
  ),
  tar_target(
    cTP.mTP_loc.exp_tp_seqs,
    read_fasta(cTP.mTP_loc.exp_tp_file)
  ),
  tar_target(
    mTP_loc.exp_tp_seqs,
    read_fasta(mTP_loc.exp_tp_file)
  ),
  tar_target(
    mTP_loc_tp.exp_seqs,
    read_fasta(mTP_loc_tp.exp_file)
  ),
  tar_target(
    SP_secreted_sp.exp_seqs,
    read_fasta(SP_secreted_sp.exp_file)
  ),
  tar_target(
    SP_secreted.er_sp.exp_seqs,
    read_fasta(SP_secreted.er_sp.exp_file)
  ),
  tar_target(
    SP_secreted.er.exp_sp_seqs,
    read_fasta(SP_secreted.er.exp_sp_file)
  ),
  tar_target(
    SP_secreted.exp_sp_seqs,
    read_fasta(SP_secreted.exp_sp_file)
  ),
  tar_target(
    SP_sp.exp_seqs,
    read_fasta(SP_sp.exp_file)
  ),
  tar_target(
    SP_sp_seqs,
    read_fasta(SP_sp_file)
  ),
  tar_target(
    TM_exp_seqs,
    read_fasta(TM_exp_file)
  ),
  tar_target(
    TM_seqs,
    read_fasta(TM_file)
  ),
  # Presequences
  tar_target(
    cTP_tp_exp_presequences,
    extract_transit_peptides(sequences = cTP_loc_tp.exp_seqs,
                             annotation_df = annotations_all,
                             remove_nonstandard = TRUE)
  ),
  tar_target(
    cTP_tp_exp,
    cTP_tp_exp_presequences[which(lengths(cTP_tp_exp_presequences) >= 10)]
  ),
  tar_target(
    cTP_tp_exp_dataset,
    write_fasta(cTP_tp_exp, paste0(data_path, "Datasets/cTP_tp_exp.fa"))
  ),
  tar_target(
    cTP_loc_exp_presequences,
    extract_transit_peptides(sequences = cTP_loc.exp_tp_seqs,
                             annotation_df = annotations_all,
                             remove_nonstandard = TRUE)
  ),
  tar_target(
    cTP_loc_exp,
    cTP_loc_exp_presequences[which(lengths(cTP_loc_exp_presequences) >= 10)]
  ),
  tar_target(
    cTP_loc_exp_dataset,
    write_fasta(cTP_loc_exp, paste0(data_path, "Datasets/cTP_loc_exp.fa"))
  ),
  tar_target(
    mTP_tp_exp_presequences,
    extract_transit_peptides(sequences = mTP_loc_tp.exp_seqs,
                             annotation_df = annotations_all,
                             remove_nonstandard = TRUE)
  ),
  tar_target(
    mTP_tp_exp,
    mTP_tp_exp_presequences[which(lengths(mTP_tp_exp_presequences) >= 10)]
  ),
  tar_target(
    mTP_tp_exp_dataset,
    write_fasta(mTP_tp_exp, paste0(data_path, "Datasets/mTP_tp_exp.fa"))
  ),
  tar_target(
    mTP_loc_exp_presequences,
    extract_transit_peptides(sequences = mTP_loc.exp_tp_seqs,
                             annotation_df = annotations_all,
                             remove_nonstandard = TRUE)
  ),
  tar_target(
    mTP_loc_exp,
    mTP_loc_exp_presequences[which(lengths(mTP_loc_exp_presequences) >= 10)]
  ),
  tar_target(
    mTP_loc_exp_dataset,
    write_fasta(mTP_loc_exp, paste0(data_path, "Datasets/mTP_loc_exp.fa"))
  ),
  tar_target(
    cTP.mTP_tp_exp_presequences,
    extract_transit_peptides(sequences = cTP.mTP_loc_tp.exp_seqs,
                             annotation_df = annotations_all,
                             remove_nonstandard = TRUE)
  ),
  tar_target(
    cTP.mTP_tp_exp,
    cTP.mTP_tp_exp_presequences[which(lengths(cTP.mTP_tp_exp_presequences) >= 10)]
  ),
  tar_target(
    cTP.mTP_tp_exp_dataset,
    write_fasta(cTP.mTP_tp_exp, paste0(data_path, "Datasets/cTP-mTP_tp_exp.fa"))
  ),
  tar_target(
    cTP.mTP_loc_exp_presequences,
    extract_transit_peptides(sequences = cTP.mTP_loc.exp_tp_seqs,
                             annotation_df = annotations_all,
                             remove_nonstandard = TRUE)
  ),
  tar_target(
    cTP.mTP_loc_exp,
    cTP.mTP_loc_exp_presequences[which(lengths(cTP.mTP_loc_exp_presequences) >= 10)]
  ),
  tar_target(
    cTP.mTP_loc_exp_dataset,
    write_fasta(cTP.mTP_loc_exp, paste0(data_path, "Datasets/cTP-mTP_loc_exp.fa"))
  ),
  tar_target(
    SP_sp_exp_presequences,
    extract_signal_peptides(sequences = SP_sp.exp_seqs, 
                            annotation_df = annotations_all,
                            remove_nonstandard = TRUE)
  ),
  tar_target(
    SP_sp_exp,
    SP_sp_exp_presequences[which(lengths(SP_sp_exp_presequences) >= 10)]
  ),
  tar_target(
    SP_sp_exp_dataset,
    write_fasta(SP_sp_exp, paste0(data_path, "Datasets/SP_sp_exp.fa"))
  ),
  tar_target(
    TM_exp_regions,
    extract_transmembrane_regions(sequences = TM_exp_seqs,
                                  annotation_df = annotations_tm,
                                  remove_nonstandard = TRUE)
  ),
  tar_target(
    TM_exp_alpha,
    TM_exp_regions[["alpha"]][which(lengths(TM_exp_regions[["alpha"]]) >= 10)]
  ),
  tar_target(
    TM_exp_beta,
    TM_exp_regions[["beta"]][which(lengths(TM_exp_regions[["beta"]]) >= 10)],
  ),
  tar_target(
    TM_exp_alpha_dataset,
    write_fasta(TM_exp_alpha, paste0(data_path, "Datasets/TM_alpha_exp.fa"))
  ),
  tar_target(
    TM_exp_beta_dataset,
    write_fasta(TM_exp_beta, paste0(data_path, "Datasets/TM_beta_exp.fa")),
  ),
  # Amyloids
  tar_target(
    amypro_file,
    paste0(data_path, "Raw_data/amypro.txt"),
    format = "file"
  ),
  tar_target(
    amypro_data,
    read_amypro_dat(amypro_file)
  ),
  tar_target(
    amypro_regions_all,
    get_amypro_regions(amypro_data)
  ),
  tar_target(
    amypro_regions,
    amypro_regions_all[which(lengths(amypro_regions_all) >= 10)]
  ),
  tar_target(
    amypro_regions_dataset,
    write_fasta(amypro_regions, paste0(data_path, "Datasets/AmyPro_regions.fa"))
  ),
  tar_target(
    cpad_file,
    paste0(data_path, "Raw_data/full_dataset_17_11_19/aggregating peptides.xlsx"),
    format = "file"
  ),
  tar_target(
    cpad_data,
    read_xlsx(cpad_file)
  ),
  tar_target(
    cpad_peptides_all,
    get_cpad_peptides(cpad_data)
  ),
  tar_target(
    cpad_peptides,
    cpad_peptides_all[which(lengths(cpad_peptides_all) >= 10)]
  ),
  tar_target(
    cpad_peptides_dataset,
    write_fasta(cpad_peptides, paste0(data_path, "Datasets/CPAD_peptides.fa"))
  ),
  tar_target(
    amyloids_combined,
    c(amypro_regions[which(!(amypro_regions %in% cpad_peptides))], cpad_peptides)
  ),
  tar_target(
    amyloids_combined_dataset,
    write_fasta(amyloids_combined, paste0(data_path, "Datasets/Amyloids_combined.fa"))
  ),
  # AMP
  tar_target(
    dbaasp_file,
    paste0(data_path, "Raw_data/dbaasp.csv"),
    format = "file"
  ),
  tar_target(
    dbaasp_data,
    read.csv(dbaasp_file)
  ),
  tar_target(
    dbaasp_amp_all,
    get_AMPs(dbaasp_data)
  ),
  tar_target(
    dbaasp_amp,
    filter_nonstandard_aa(dbaasp_amp_all[which(lengths(dbaasp_amp_all) >= 10)])
  ),
  tar_target(
    dbaasp_amp_max100,
    dbaasp_amp[which(lengths(dbaasp_amp) <= 100)]
  ),
  tar_target(
    dbaasp_amp_dataset,
    write_fasta(dbaasp_amp, paste0(data_path, "Datasets/DBAASP_AMP.fa"))
  ),
  tar_target(
    dbaasp_amp_max100_dataset,
    write_fasta(dbaasp_amp_max100, paste0(data_path, "Datasets/DBAASP_AMP_max100.fa"))
  ),
  tar_target(
    amyloids_combined_70,
    filter_with_cdhit(amyloids_combined, 0.7)
  ),
  tar_target(
    amypro_regions_70,
    filter_with_cdhit(amypro_regions, 0.7)
  ),
  tar_target(
    cpad_peptides_70,
    filter_with_cdhit(cpad_peptides, 0.7)
  ),
  tar_target(
    cTP_loc_exp_70,
    filter_with_cdhit(cTP_loc_exp, 0.7)
  ),
  tar_target(
    cTP.mTP_loc_exp_70,
    filter_with_cdhit(cTP.mTP_loc_exp, 0.7)
  ),
  tar_target(
    cTP.mTP_tp_exp_70,
    filter_with_cdhit(cTP.mTP_tp_exp, 0.7)
  ),
  tar_target(
    cTP_tp_exp_70,
    filter_with_cdhit(cTP_tp_exp, 0.7)
  ),
  tar_target(
    dbaasp_amp_70,
    filter_with_cdhit(dbaasp_amp, 0.7)
  ),
  tar_target(
    dbaasp_amp_max100_70,
    filter_with_cdhit(dbaasp_amp_max100, 0.7)
  ),
  tar_target(
    mTP_loc_exp_70,
    filter_with_cdhit(mTP_loc_exp, 0.7)
  ),
  tar_target(
    mTP_tp_exp_70,
    filter_with_cdhit(mTP_tp_exp, 0.7)
  ),
  tar_target(
    SP_sp_exp_70,
    filter_with_cdhit(SP_sp_exp, 0.7)
  ),
  tar_target(
    TM_exp_alpha_70,
    filter_with_cdhit(TM_exp_alpha, 0.7)
  ),
  tar_target(
    TM_exp_beta_70,
    filter_with_cdhit(TM_exp_beta, 0.7)
  ),
  tar_target(
    datasets_list,
    list("Amyloids combined" = amyloids_combined, "AmyPro regions" = amypro_regions, "CPAD peptides" = cpad_peptides,
         "cTP experimentally verified location" = cTP_loc_exp, "cTP-mTP experimentally verified location" = cTP.mTP_loc_exp,
         "cTP-mTP experimentally verified presequence" = cTP.mTP_tp_exp, "cTP experimentally verified presequence" = cTP_tp_exp,
         "DBAASP AMP" = dbaasp_amp, "DBAASP AMP max 100 aa" = dbaasp_amp_max100, "mTP experimentally verified location" = mTP_loc_exp,
         "mTP experimentally verified presequence" = mTP_tp_exp, "SP experimentally verified presequence" = SP_sp_exp,
         "TM regions experimentally verified - alpha" = TM_exp_alpha, "TM regions experimentally verified - beta" = TM_exp_beta)
  ),
  tar_target(
    datasets_list_cdhit,
    list("Amyloids combined" = amyloids_combined_70, "AmyPro regions" = amypro_regions_70, "CPAD peptides" = cpad_peptides_70,
         "cTP experimentally verified location" = cTP_loc_exp_70, "cTP-mTP experimentally verified location" = cTP.mTP_loc_exp_70,
         "cTP-mTP experimentally verified presequence" = cTP.mTP_tp_exp_70, "cTP experimentally verified presequence" = cTP_tp_exp_70,
         "DBAASP AMP" = dbaasp_amp_70, "DBAASP AMP max 100 aa" = dbaasp_amp_max100_70, "mTP experimentally verified location" = mTP_loc_exp_70,
         "mTP experimentally verified presequence" = mTP_tp_exp_70, "SP experimentally verified presequence" = SP_sp_exp_70,
         "TM regions experimentally verified - alpha" = TM_exp_alpha_70, "TM regions experimentally verified - beta" = TM_exp_beta_70)
  )
)

