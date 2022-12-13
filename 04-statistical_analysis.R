library(dplyr)
library(biogram)
library(seqR)
library(targets)
library(purrr)

source("./functions/extract_features.R")
source("./functions/motif_functions.R")

tar_load(c(datasets_list, taxonomic_representation))
df_freq_tax2 <- readRDS("df_freq_tax2.rds") 

data_path <- "/media/kasia/Data/Dropbox/Presequences/"

filtered_datasets_list <- datasets_list[c("cTP experimentally verified presequence", 
                                          "mTP experimentally verified presequence", 
                                          "SP experimentally verified presequence")]


k_vec <- c(rep(2, 4), rep(3, 16), rep(4, 27), rep(5, 16))
gap_vec <- c(0:3, c(0, 0), c(0, 1), c(1, 0), c(1, 1), 
             c(0, 2), c(2, 0), c(1, 2), c(2, 1), c(2, 2), 
             c(0, 3), c(3, 0), c(1, 3), c(3, 1), c(2, 3), 
             c(3, 2), c(3, 3),
             c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), 
             c(1, 1, 0), c(1, 0, 1), c(0, 1, 1), c(1, 1, 1), 
             c(2, 0, 0), c(0, 2, 0), c(0, 0, 2), c(2, 2, 0), 
             c(2, 0, 2), c(0, 2, 2), c(2, 2, 2), c(2, 1, 0), 
             c(2, 0, 1), c(1, 2, 0), c(1, 0, 2), c(0, 1, 2), 
             c(0, 2, 1), c(2, 1, 1), c(1, 2, 1), c(1, 1, 2), 
             c(1, 2, 2), c(2, 1, 2), c(2, 2, 1),
             c(0, 0, 0, 0), c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 1, 0), 
             c(0, 0, 0, 1), c(1, 1, 0, 0), c(1, 0, 1, 0), c(1, 0, 0, 1), 
             c(0, 1, 1, 0), c(0, 1, 0, 1), c(0, 0, 1, 1), c(1, 1, 1, 0), 
             c(1, 0, 1, 1), c(1, 1, 0, 1), c(0, 1, 1, 1), c(1, 1, 1, 1))


motif_data <- df_freq_tax2 %>% 
  select(c(Type, Dataset, `Frequent in`, Motif)) %>% 
  mutate(Motif = gsub(" ", "", Motif)) %>% 
  unique()


ngram_data <- mapply(function(k, gap) {
  lapply(names(filtered_datasets_list), function(ith_set) {
    bin_matrix <- calculate_binary_ngram_matrix(filtered_datasets_list[[ith_set]], k, gap)
    colnames(bin_matrix) <- decode_ngrams(colnames(bin_matrix))
    ngram_data <- bin_matrix[, which(colnames(bin_matrix) %in% motif_data[["Motif"]])] %>% 
      mutate(dataset = strsplit(ith_set, " ")[[1]][1],
             Entry = names(filtered_datasets_list[[ith_set]])) %>% 
      left_join(select(taxonomic_representation, c("Entry", "superkingdom", "kingdom", "phylum", "class")))
  }) %>% reduce(full_join)
}, k_vec[1:5], gap_vec[1:5]) %>% reduce(full_join) %>% 
  mutate(kingdom = ifelse(is.na(kingdom), "Unknown", kingdom),
         phylum = ifelse(is.na(phylum), "Unknown", phylum),
         class = ifelse(is.na(class), "Unknown", class))

test_res <- lapply(unique(motif_data[["Dataset"]]), function(ith_dataset) {
  lapply(unique(filter(motif_data, Dataset == ith_dataset)[["Frequent in"]]), function(ith_frequent) {
    x <- ngram_data %>%
      filter(dataset == strsplit(ith_dataset, " ")[[1]][1]) %>% 
      .[, which(colnames(ngram_data) %in% c(unname(filter(motif_data, Dataset == ith_dataset, `Frequent in` == ith_frequent)[["Motif"]]),
                                            "Entry", "superkingdom", "kingdom", "phylum", "class"))]
    tax_group <- case_when(ith_frequent %in% c("Eukaryota", "Bacteria", "Archaea", "Viruses") ~ "superkingdom",
                           ith_frequent %in% c("Viridiplantae", "Fungi", "Metazoa") ~ "kingdom",
                           ith_frequent %in% c("Streptophyta", "Chlorophyta", "Chordata", "Arthropoda") ~ "phylum",
                           ith_frequent %in% c("Mammalia", "Insecta", "Lepidosauria", "Arachnida") ~ "class")
    groups <- if(tax_group == "phylum" & ith_dataset == "cTP") {
      c("Streptophyta", "Chlorophyta")
    } else if(ith_dataset == "SP (phylum)") {
      c("Chordata", "Arthropoda")
    } else if(tax_group == "kingdom" & ith_dataset == "cTP") {
      c("Viridiplantae", "Unknown")
    } else if(tax_group == "kingdom" & ith_dataset != "cTP") {
      c("Viridiplantae", "Fungi", "Metazoa")
    } else if(tax_group == "superkingdom") {
      c("Eukaryota", "Bacteria", "Archaea", "Viruses") 
    } else if(tax_group == "class") {
      c("Mammalia", "Insecta", "Lepidosauria", "Arachnida")
    }
    print(paste0(ith_dataset))
    print(paste0(ith_frequent))
    print(paste0(tax_group))
    print(paste0(groups))
    combns <- combn(groups, 2, simplify = FALSE)
    lapply(1:length(combns), function(i) {
      lapply(colnames(x[which(!(colnames(x) %in% c("Entry", "superkingdom", "kingdom", "phylum", "class")))]), function(ith_motif) {
        print(paste0("Dataset: ", ith_dataset, " --- Frequent in: ", ith_frequent, " --- Motif: ", ith_motif))
        data.frame(
          Dataset = ith_dataset,
          Frequent = ith_frequent,
          Motif = ith_motif,
          Tax1 = combns[[i]][1],
          Tax2 = combns[[i]][2],
          pval = wilcox.test(filter(x, get(tax_group) == combns[[i]][1])[[ith_motif]],
                             filter(x, get(tax_group) == combns[[i]][2])[[ith_motif]])[["p.value"]]
        )
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
}) %>% bind_rows() 



