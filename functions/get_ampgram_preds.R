get_ampgram_preds <- function(data_path, organism, nthreads, seq_type) {
  seqs <- read_fasta(paste0(data_path, "Data/Sequences_for_analysis/", name, "_", seq_type, ".fa"))
  preds <- pblapply(seq(1, length(seqs), by = 500), cl = nthreads, function(i) {
    if(i + 499 > length(seqs)) {
      pred2df(predict(AmpGram_model, seqs[i:length(seqs)]))
    } else {
      pred2df(predict(AmpGram_model, seqs[i:(i+499)]))
    }
  }) %>% bind_rows()
  write.csv(preds, paste0(data_path, "Analizy/Wyniki/", organism, "_", seq_type, "_AmpGram_results.csv"),
            row.names = FALSE)
}
