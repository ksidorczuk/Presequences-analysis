# This script uses datasets obtained in a targets pipeline 
# and runs prediction of their antimicrobial properties
library(dplyr)
library(biogram)
library(targets)
library(ampir)
library(AmpGramModel)
library(seqR)
library(ranger)

tar_load(list.files("_targets/objects"))
source("functions/predict_AMPs.R")
source("functions/predict_with_AmpGram.R")

data_path <- "~/Dropbox/Presequences/Datasets"
seq_files <- list.files(data_path, full.names = TRUE)


lapply(seq_files, function(ith_file) {
  out_path <- gsub(".fa", "", gsub("Datasets", "Prediction_results", ith_file))
  print(paste0("Running ampir predictions for file: ", ith_file))
  predict_with_ampir(ith_file, paste0(out_path, "_ampir.csv"))
  print(paste0("Running amPEPpy predictions for file: ", ith_file))
  predict_with_amPEPpy(ith_file, paste0(out_path, "_amPEPpy.tsv"))
  print(paste0("Running AmpGram predictions for file: ", ith_file))
  predict_with_AmpGram(ith_file, gsub("Datasets", "Prediction_results/", data_path), paste0(out_path, "_AmpGram.csv"))
})

lapply(seq_files, function(ith_file) {
  predict_with_ampscanner(ith_file, "~/AMP_Scanner2_Feb2020_Model/AMP_Scanner2_Feb2020_Model/")
})