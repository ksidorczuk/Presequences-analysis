process_ampeppy_res_files <- function(res_file) {
  df <- read.delim(res_file) %>% 
    select(-probability_nonAMP) %>% 
    setNames(c("Probability", "Prediction", "Seq_ID"))
  x <- last(strsplit(res_file, "/")[[1]])
  mutate(df,
         Probability = as.numeric(Probability),
         Seq_ID = as.character(df[["Seq_ID"]]),
         Dataset = gsub("_amPEPpy.tsv", "", x),
         Software = "amPEPpy",
         Prediction = ifelse(Prediction == "AMP", TRUE, FALSE))
}

process_ampgram_res_files <- function(res_file) {
  df <- read.csv(res_file) %>% 
    setNames(c("Seq_ID", "Probability", "Prediction"))
  x <- last(strsplit(res_file, "/")[[1]])
  mutate(df,
         Probability = as.numeric(Probability),
         Seq_ID = as.character(df[["Seq_ID"]]),
         Dataset = gsub("_AmpGram.csv", "", x),
         Software = "AmpGram")
}

process_ampscanner_res_files <- function(res_file) {
  df <- read.csv(res_file) %>% 
    select(-Sequence) %>% 
    setNames(c("Seq_ID", "Prediction", "Probability"))
  x <- last(strsplit(res_file, "/")[[1]])
  mutate(df,
         Probability = as.numeric(Probability),
         Seq_ID = as.character(df[["Seq_ID"]]),
         Dataset = gsub("_AMPScanner.csv", "", x),
         Software = "AmpScanner",
         Prediction = ifelse(Prediction == "AMP", TRUE, FALSE))
}

process_ampir_res_files <- function(res_file) {
  df <- read.csv(res_file) %>% 
    select(-seq_aa) %>% 
    setNames(c("Seq_ID", "Probability")) %>% 
    mutate(Prediction = ifelse(Probability > 0.5, TRUE, FALSE))
  x <- last(strsplit(res_file, "/")[[1]])
  mutate(df,
         Probability = as.numeric(Probability),
         Seq_ID = as.character(df[["Seq_ID"]]),
         Dataset = gsub("_ampir.csv", "", x),
         Software = "ampir") 
}

