library(dplyr)
library(biogram)
library(seqR)
library(targets)
library(purrr)
library(ggplot2)
library(tidyr)
library(ggdendro)

source("./functions/extract_features.R")
source("./functions/motif_functions.R")

tar_load(c(datasets_list, taxonomic_representation))
df_freq_tax2 <- readRDS("df_freq_tax2.rds") 

data_path <- "/media/kasia/Data/Dropbox/Presequences/"

filtered_datasets_list <- datasets_list[c("cTP experimentally verified presequence", 
                                          "mTP experimentally verified presequence", 
                                          "SP experimentally verified presequence")]


k_vec <- c(rep(list(2), 4), rep(list(3), 16), rep(list(4), 27), rep(list(5), 16))
gap_vec <- list(0, 1, 2, 3, list(c(0, 0)), list(c(0, 1)), list(c(1, 0)), list(c(1, 1)), 
                list(c(0, 2)), list(c(2, 0)), list(c(1, 2)), list(c(2, 1)), list(c(2, 2)), 
                list(c(0, 3)), list(c(3, 0)), list(c(1, 3)), list(c(3, 1)), list(c(2, 3)), 
                list(c(3, 2)), list(c(3, 3)),
                list(c(0, 0, 0)), list(c(1, 0, 0)), list(c(0, 1, 0)), list(c(0, 0, 1)), 
                list(c(1, 1, 0)), list(c(1, 0, 1)), list(c(0, 1, 1)), list(c(1, 1, 1)), 
                list(c(2, 0, 0)), list(c(0, 2, 0)), list(c(0, 0, 2)), list(c(2, 2, 0)), 
                list(c(2, 0, 2)), list(c(0, 2, 2)), list(c(2, 2, 2)), list(c(2, 1, 0)), 
                list(c(2, 0, 1)), list(c(1, 2, 0)), list(c(1, 0, 2)), list(c(0, 1, 2)), 
                list(c(0, 2, 1)), list(c(2, 1, 1)), list(c(1, 2, 1)), list(c(1, 1, 2)), 
                list(c(1, 2, 2)), list(c(2, 1, 2)), list(c(2, 2, 1)),
                list(c(0, 0, 0, 0)), list(c(1, 0, 0, 0)), list(c(0, 1, 0, 0)), list(c(0, 0, 1, 0)), 
                list(c(0, 0, 0, 1)), list(c(1, 1, 0, 0)), list(c(1, 0, 1, 0)), list(c(1, 0, 0, 1)), 
                list(c(0, 1, 1, 0)), list(c(0, 1, 0, 1)), list(c(0, 0, 1, 1)), list(c(1, 1, 1, 0)), 
                list(c(1, 0, 1, 1)), list(c(1, 1, 0, 1)), list(c(0, 1, 1, 1)), list(c(1, 1, 1, 1)))


motif_data <- df_freq_tax2 %>% 
  select(c(Type, Dataset, `Frequent in`, Motif)) %>% 
  mutate(Motif = gsub(" ", "", Motif)) %>% 
  unique()


ngram_data <- mapply(function(k, gap) {
  print(paste0("k: ", k))
  print(paste0("gap: ", gap))
  lapply(names(filtered_datasets_list), function(ith_set) {
    bin_matrix <- calculate_binary_ngram_matrix(filtered_datasets_list[[ith_set]], k, gap)
    colnames(bin_matrix) <- decode_ngrams(colnames(bin_matrix))
    to_select <- colnames(bin_matrix)[which(colnames(bin_matrix) %in% motif_data[["Motif"]])]
    ngram_dat <- if(length(to_select) >= 1) {
      select(bin_matrix, to_select) %>%
        mutate(dataset = strsplit(ith_set, " ")[[1]][1],
               Entry = names(filtered_datasets_list[[ith_set]])) %>%
        left_join(., select(taxonomic_representation, c("Entry", "superkingdom", "kingdom", "phylum", "class")))
    } else {
      data.frame(dataset = strsplit(ith_set, " ")[[1]][1],
                 Entry = names(filtered_datasets_list[[ith_set]])) %>%
        left_join(., select(taxonomic_representation, c("Entry", "superkingdom", "kingdom", "phylum", "class")))
    }
    print(paste0("Dataset: ", ith_set, " --- k: ", k, " --- gaps: ", gap))
    ngram_dat
  }) %>% reduce(full_join)
}, k_vec, gap_vec, SIMPLIFY = FALSE) %>% reduce(full_join) %>% 
  mutate(kingdom = ifelse(is.na(kingdom), "Unknown", kingdom),
         phylum = ifelse(is.na(phylum), "Unknown", phylum),
         class = ifelse(is.na(class), "Unknown", class))

ngram_data[is.na(ngram_data)] <- 0
saveRDS(ngram_data, "ngram_data.rds")

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
      y <- filter(x, get(tax_group) %in% combns[[i]])
      y <- cbind(select(y, c("Entry", "superkingdom", "kingdom", "phylum", "class")), 
                 select(y, !c("Entry", "superkingdom", "kingdom", "phylum", "class"))[, which(colSums(select(y, !c("Entry", "superkingdom", "kingdom", "phylum", "class"))) > 0)])
      lapply(colnames(y[which(!(colnames(y) %in% c("Entry", "superkingdom", "kingdom", "phylum", "class")))]), function(ith_motif) {
        print(paste0("Dataset: ", ith_dataset, " --- Frequent in: ", ith_frequent, " --- Motif: ", ith_motif))
        data.frame(
          Dataset = ith_dataset,
          # Frequent = ith_frequent,
          Motif = ith_motif,
          Tax1 = combns[[i]][1],
          Tax2 = combns[[i]][2],
          pval = wilcox.test(filter(y, get(tax_group) == combns[[i]][1])[[ith_motif]],
                             filter(y, get(tax_group) == combns[[i]][2])[[ith_motif]])[["p.value"]]
        )
      }) %>% bind_rows() 
    }) %>% bind_rows()
  }) %>% bind_rows()
}) %>% bind_rows() %>% 
  unique()
#4211
test_res %>% 
  mutate(Comparison = paste0(Tax1, ' vs. ', Tax2),
         is_significant = ifelse(pval < 0.05, TRUE, FALSE)) %>% 
  ggplot(aes(y = Comparison, x = Motif, fill = is_significant)) +
  geom_tile() +
  facet_wrap(~Dataset, scales = "free", ncol = 1) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90))


plot_dat <- test_res %>% 
  mutate(Comparison = paste0(Tax1, ' vs. ', Tax2),
         is_significant = ifelse(pval < 0.05, TRUE, FALSE)) %>% 
  select(c("Dataset", "Motif", "Comparison", "is_significant"))

plot_dat_wide <- pivot_wider(plot_dat, c(Dataset, Comparison), names_from = Motif, values_from = is_significant)

lapply(unique(plot_dat[["Dataset"]]), function(ith_set) {
  p <- plot_dat %>% 
    filter(Dataset == ith_set) %>% 
    ggplot(aes(y = Comparison, x = Motif, fill = is_significant)) +
    geom_tile(color = "white") +
    ggtitle(ith_set) +
    theme_bw(base_size = 5) +
    coord_fixed(ratio = 1) +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "bottom")
  ggsave(paste0(data_path, "ngram_results/Statistical_analysis_", gsub(" ", "_", ith_set), ".png"),
         p, width = 6+nrow(filter(plot_dat, Dataset == ith_set))/70, height = 4, limitsize = FALSE)
})


