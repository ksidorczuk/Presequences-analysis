# This scripts contain code used for analyses of proteins
# imported via endomembrane system.

library(biogram)
library(dplyr)
library(ggplot2)

source("functions/aa_comp_functions.R")
source("functions/extract_features.R")

enc10 <- list("s" = "S", 
              "t" = c("T", "N", "Q"), 
              "a" = "A",
              "v" = c("V", "I", "L"),
              "m" = "M",
              "f" = c("F", "Y", "W"),
              "g" = c("C", "G", "P"),
              "d" = c("D", "E"),
              "r" = "R",
              "h" = c("H", "K"))

enc7 <- list("s" = c("S", "T", "N", "Q"), 
             "v" = c("A", "V", "I", "L"),
             "m" = "M",
             "f" = c("F", "Y", "W"),
             "g" = c("C", "G", "P"),
             "d" = c("D", "E"),
             "r" = c("R", "H", "K"))

enc8 <- list("s" = c("S", "T"),
             "n" = c("N", "Q"), 
             "v" = c("A", "V", "I", "L"),
             "m" = "M",
             "f" = c("F", "Y", "W"),
             "g" = c("C", "G", "P"),
             "d" = c("D", "E"),
             "r" = c("R", "H", "K"))

SDM12 <- list("a" =	"A",
              "d" = "D",
              "k" = c("K", "E", "R"),
              "n" = "N",
              "t" = c("T", "S", "Q"),
              "y" = c("Y", "F"),
              "l" = c("L", "I", "V", "M"), 
              "c" = "C",
              "w" = "W",
              "h" = "H",
              "g" = "G",
              "p" = "P")

datasets <- c("AMY", "AMY_withoutAMYC1", "CAH", "NPP", "SOD")
alphabets <- c("enc10", "enc8", "enc7", "SDM12")


# Translate sequences using reduced alphabets
lapply(datasets, function(ith_set) {
  lapply(alphabets, function(ith_alphabet) {
    s <- read_fasta(paste0("/media/kasia/Data/Dropbox/Presequences/Proteins improted via endomembrane system/comparisons/", ith_set, ".fa"))
    rs <- sapply(s, function(i) degenerate(i, get(ith_alphabet)))
    write_fasta(rs, paste0("/media/kasia/Data/Dropbox/Presequences/Proteins improted via endomembrane system/comparisons/", ith_set, "_", ith_alphabet, ".fa"))
  })
})


# Read in original sequence datasets for comparisons
all_datasets <- lapply(datasets, function(ith_set) {
  read_fasta(paste0("/media/kasia/Data/Dropbox/Presequences/Proteins improted via endomembrane system/comparisons/", ith_set, ".fa")) 
}) %>% setNames(datasets)

aa_comp <- calculate_aa_comp_peptides(all_datasets) %>% 
  mutate(`Imported via ES` = factor(ifelse(prot %in% c("AMY1", "AMY3D", "ATCA1", "NPP1", "NPP2", "NPP6", "sp|Q43008|SODM"), TRUE, FALSE)),
         dataset = ifelse(dataset == "AMY_withoutAMYC1", "AMY (without AMYC1)", dataset))

# Compare amino acid composition
lapply(unique(aa_comp[["dataset"]]), function(ith_set) {
  p <- filter(aa_comp, dataset == ith_set) %>% 
    ggplot(aes(x = prot, y = Frequency, fill = `Imported via ES`, group = prot)) +
    facet_wrap(~`Amino acid`) +
    geom_col(position = "dodge") +
    ggtitle(ith_set) +
    coord_flip() +
    theme_bw() +
    scale_fill_manual("Imported via endomembrane system", values = c("TRUE" = "#e49144", "FALSE" = "grey70")) +
    theme(legend.position = "bottom") +
    xlab("Protein")
  ggsave(plot = p, filename = paste0("/media/kasia/Data/Dropbox/Presequences/Proteins improted via endomembrane system/comparisons/", ith_set, "_aa_comp.eps"),
           width = 12, height = 9)
})

# Compare charge distribution
charge_dat <- lapply(names(all_datasets), function(ith_set){
  x <- all_datasets[[ith_set]]
  lapply(names(all_datasets[[ith_set]]), function(i){
    data.frame(Dataset = ith_set,
               Protein = i, 
               aa = x[[i]], 
               pos = 1:length(x[[i]])) %>% 
      mutate(charge = factor(case_when(aa %in% c("K", "R", "H") ~ "Positive",
                                aa %in% c("D", "E") ~ "Negative",
                                aa %in% c("A", "C", "F", "G", "I", "L", "M", "N", 
                                          "P", "Q", "S", "T", "V", "W", "Y") ~ "Neutral"), 
                             levels = c("Negative", "Neutral", "Positive")))
  }) %>% bind_rows()
}) %>% bind_rows() %>% 
  mutate(`Imported via ES` = factor(ifelse(Protein %in% c("AMY1", "AMY3D", "ATCA1", "NPP1", "NPP2", "NPP6", "sp|Q43008|SODM"), TRUE, FALSE)),
         dataset = ifelse(Dataset == "AMY_withoutAMYC1", "AMY (without AMYC1)", Dataset))

lapply(unique(charge_dat[["Dataset"]]), function(ith_set) {
  p <- filter(charge_dat, Dataset == ith_set) %>% 
    ggplot(aes(x = pos, y = charge, color = `Imported via ES`, group = Protein)) +
    geom_line() +
    facet_wrap(~Protein, ncol = 1) +
    ggtitle(ith_set) +
    theme_bw() +
    xlab("position") +
    scale_color_manual("Imported via endomembrane system", values = c("TRUE" = "#e49144", "FALSE" = "grey70")) +
    theme(legend.position = "bottom")
  ggsave(plot = p, filename = paste0("/media/kasia/Data/Dropbox/Presequences/Proteins improted via endomembrane system/comparisons/", ith_set, "_charge.eps"),
         width = 15, height = 8)
})
