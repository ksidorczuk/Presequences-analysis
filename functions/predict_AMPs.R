#' Run AMP predictions using amPEPpy
#' 
#' This function reads sequences of peptides from a file and runs
#' predictions using amPEPpy software. 
#' @param sequence_file name of a file containing sequences in a FASTA
#' format to analyse
#' @param outfile name, and optionally path, of the output file
#' @param n_threads number of threads used for running amPEPpy, 
#' by default 8 threads are used
#' @param amPEPpy_model_path path to the amPEP.model file required
#' for running predictions
#' @return NULL
predict_with_amPEPpy <- function(sequence_file, out_file, n_threads = 8, amPEPpy_model_path = "~/Programy/amPEPpy/amPEP.model") {
  system(paste0("ampep predict -m ", amPEPpy_model_path, " -i ", sequence_file, " -o ", out_file, " -t ", n_threads))
}

#' Run AMP predictions using AMPlify
#' 
#' This function reads sequences of peptides from a file and runs
#' predictions using AMPlify software.
#' @param sequence_file name of a file containing sequences in a FASTA
#' format to analyse
#' @param outfile name, and optionally path, of the output file
#' @param AMPlify_path path to the AMPlify main directory
#' @return NULL
predict_with_AMPlify <- function(sequence_file, outfile, AMPlify_path = "~/Programy/AMPlify-1.0.0/") {
  system(paste0("python3 ", AMPlify_path, "src/AMPlify.py -md ", AMPlify_path, "models/ -s ", sequence_file, " -of 'tsv'"))
  file.rename(list.files(pattern = "AMPlify"), outfile)
}

#' Run AMP predictions using ampir
#' 
#' This function reads sequences of peptides from a file and runs
#' predictions using ampir package. It uses chunks of 20,000 sequences
#' at once to allow analysis of large files step by step. All results 
#' are then combined.
#' @param sequence_file name of a file containing sequences in a FASTA
#' format to analyse
#' @param outfile name, and optionally path, of the output file
#' @return NULL
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

#' Run AMP predictions using AMPScanner
#' 
#' This function reads sequences of peptides from a file and runs
#' predictions using AMPScanner V2 software.
#' @param sequence_file name of a file containing sequences in a FASTA
#' format to analyse
#' @param ampscanner_path path to the AMPScanner main directory
#' @return NULL
predict_with_ampscanner <- function(sequence_file, ampscanner_path) {
  system(paste0("python3 ", ampscanner_path, "amp_scanner_predict.py ", sequence_file, " ", ampscanner_path, "TrainedModels/021820_FULL_MODEL.h5"))
}