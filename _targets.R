library(targets)
library(dplyr)
library(readxl)
library(biogram)

data_path <- "~/Downloads/Datasets/"

source("functions/extract_presequences.R")
source("functions/extract_transmembrane_regions.R")
source("functions/get_sequence_numbers.R")

list(
  # Annotation files
  tar_target(
    annotations_all_file,
    paste0(data_path, "uniprot-reviewed yes.xlsx"),
    format = "file"
  ),
  tar_target(
    annotations_tm_file,
    paste0(data_path, "uniprot-annotation (type transmem)-filtered-reviewed yes.xlsx"),
    format = "file"
  ),
  tar_target(
    annotations_intm_file,
    paste0(data_path, "uniprot-annotation (type intramem)+reviewed yes.xlsx"),
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
    paste0(data_path, "cTP_loc_tp-exp.fasta"),
    format = "file"
  ),
  tar_target(
    cTP_loc.exp_tp_file,
    paste0(data_path, "cTP_loc-exp_tp.fasta"),
    format = "file"
  ),
  tar_target(
    `cTP.mTP_loc_exp.tp_file`,
    paste0(data_path, "cTP-mTP_loc_exp-tp.fasta"),
    format = "file"
  ),
  tar_target(
    `cTP.mTP_loc.exp_tp_file`,
    paste0(data_path, "cTP-mTP_loc-exp_tp.fasta"),
    format = "file"
  ),
  tar_target(
    `mTP_loc.exp_tp_file`,
    paste0(data_path, "mTP_loc-exp_tp.fasta"),
    format = "file"
  ),
  tar_target(
    `mTP_loc_tp.exp_file`,
    paste0(data_path, "mTP_loc_tp-exp.fasta"),
    format = "file"
  ),
  tar_target(
    `SP_secreted_sp.exp_file`,
    paste0(data_path, "SP_secreted_sp-exp.fasta"),
    format = "file"
  ),
  tar_target(
    `SP_secreted.er_sp.exp_file`,
    paste0(data_path, "SP_secreted-er_sp-exp.fasta"),
    format = "file"
  ),
  tar_target(
    `SP_secreted.er.exp_sp_file`,
    paste0(data_path, "SP_secreted-er-exp_sp.fasta"),
    format = "file"
  ),
  tar_target(
    `SP_secreted.exp_sp_file`,
    paste0(data_path, "SP_secreted-exp_sp.fasta"),
    format = "file"
  ),
  tar_target(
    `SP_sp.exp_file`,
    paste0(data_path, "SP_sp-exp.fasta"),
    format = "file"
  ),
  tar_target(
    `SP_sp_file`,
    paste0(data_path, "SP_sp.fasta"),
    format = "file"
  ),
  tar_target(
    `TM_exp_file`,
    paste0(data_path, "TM_exp.fasta"),
    format = "file"
  ),
  tar_target(
    `TM_file`,
    paste0(data_path, "TM.fasta"),
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
    cTP.mTP_loc_exp.tp_seqs,
    read_fasta(cTP.mTP_loc_exp.tp_file)
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
  )
)

