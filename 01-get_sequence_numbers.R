library(targets)
library(dplyr)
library(biogram)

tar_load(list.files("_targets/objects"))

TP_summ <- summarise_TP_info(c("cTP_loc_tp.exp_seqs", "cTP_loc.exp_tp_seqs", "cTP.mTP_loc_exp.tp_seqs", 
                               "cTP.mTP_loc.exp_tp_seqs", "mTP_loc_tp.exp_seqs", "mTP_loc.exp_tp_seqs"),
                             annotations_all)

SP_summ <- summarise_SP_info(c("SP_secreted_sp.exp_seqs", "SP_secreted.er_sp.exp_seqs", "SP_secreted.er.exp_sp_seqs", 
                               "SP_secreted.exp_sp_seqs", "SP_sp_seqs", "SP_sp.exp_seqs"),
                             annotations_all)

TM_summ <- summarise_TP_info(c("TM_exp_seqs", "TM_seqs"),
                             annotations_tm)
