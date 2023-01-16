#' Filter with CD-HIT
#' 
#' This functions uses external software CD-HIT to cluster sequences 
#' into clusters that meet a given sequence identity threshold. Each 
#' cluster has one representative sequence. This function returns 
#' a list of sequences that represent each cluster.
#' 
#' @param input_seq list of input sequences
#' @param thresh \code{double} indicating sequence identity threshold 
#' (e.g. 0.7 means 70% identity)
#' @param word_length word size
#' @param cdhit_path \code{character} path to cd-hit
#' @return list of filtered sequences
filter_with_cdhit <- function(sequences, threshold, word_length = 2, cdhit_path = "./third-party") {
  input <- tempfile(tmpdir = getwd())
  output <- tempfile(tmpdir = getwd())
  cdhit <- paste0(cdhit_path, "/cdhit -i ", input,  " -o ", output, " -c ", threshold, " -n ", word_length)
  write_fasta(sequences, input)
  system(cdhit)
  res <- read_fasta(output)
  file.remove(input, output, paste0(output, ".clstr"))
  res
}
