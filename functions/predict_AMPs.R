predict_with_amPEPpy <- function(sequence_file, out_file, n_threads = 8, amPEPpy_model_path = "~/Programy/amPEPpy/amPEP.model") {
  system(paste0("ampep predict -m ", amPEPpy_model_path, " -i ", sequence_file, " -o ", out_file, " -t ", n_threads))
}

predict_with_AMPlify <- function(sequence_file, outfile, AMPlify_path = "~/Programy/AMPlify-1.0.0/") {
  system(paste0("python3 ", AMPlify_path, "src/AMPlify.py -md ", AMPlify_path, "models/ -s ", sequence_file, " -of 'tsv'"))
  file.rename(list.files(pattern = "AMPlify"), outfile)
}

predict_with_ampir <- function(sequence_file, outfile) {
  seqs <- read_faa(sequence_file) 
  
  preds <- lapply(seq(1, nrow(seqs), 20000), function(i) {
    x <- if(i+19999 < nrow(seqs)) {
      i+19999
    } else {
      nrow(seqs)
    }
    predict_amps(seqs[i:x, ], model = "mature", n_cores = 8)
  }) %>% bind_rows()
  write.csv(preds, file = paste0(outfile),
            row.names = FALSE)
}

predict_with_ampscanner <- function(sequence_file, ampscanner_path) {
  system(paste0("python3 ", ampscanner_path, "amp_scanner_predict.py ", sequence_file, " ", ampscanner_path, "TrainedModels/021820_FULL_MODEL.h5"))
}