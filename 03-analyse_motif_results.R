library(dplyr)
library(ggplot2)
library(readxl)
library(ggrepel)
library(ggbeeswarm)
library(biogram)
library(tidyr)
library(tidytext)
library(ggtext)
library(patchwork)

data_path <- "/media/kasia/Data/Dropbox/Presequences/"

decode_motifs <- function(split_motif, encoding) {
  if(encoding == "enc7") {
    plyr::mapvalues(split_motif, from=c("1", "3", "5", "6", "7", "8", "9"), to=c("[STNQ]", "[AVIL]", "M", "[FYW]", "[CGP]", "[DE]", "[RHK]"))
  } else if(encoding == "enc8") {
    plyr::mapvalues(split_motif, from=c("1", "2", "3", "5", "6", "7", "8", "9"), to=c("[ST]", "[NQ]", "[AVIL]", "M", "[FYW]", "[CGP]", "[DE]", "[RHK]"))
  } else if(encoding == "enc10") {
    plyr::mapvalues(split_motif, from=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j"), to=c("S", "[TNQ]", "A", "[VIL]", "M", "[FYW]", "[CGP]", "[DE]", "R", "[HK]"))
  }
}

color_encodings <- function(x) {
  x %>% 
    mutate(Motif = gsub("[AVIL]", "<span style='color:#4a219c'>[AVIL]</span>", Motif, fixed = TRUE)) %>% 
    mutate(Motif = gsub("[VIL]", "<span style='color:#4a219c'>[VIL]</span>", Motif, fixed = TRUE)) %>% 
    mutate(Motif = gsub("[CGP]", "<span style='color:#9c9021'>[CGP]</span>", Motif, fixed = TRUE)) %>% 
    mutate(Motif = gsub("[STNQ]", "<span style='color:#21739c'>[STNQ]</span>", Motif, fixed = TRUE)) %>% 
    mutate(Motif = gsub("[TNQ]", "<span style='color:#21739c'>[TNQ]</span>", Motif, fixed = TRUE)) %>% 
    mutate(Motif = gsub("[ST]", "<span style='color:#21739c'>[ST]</span>", Motif, fixed = TRUE)) %>% 
    mutate(Motif = gsub("[FYW]", "<span style='color:#8e219c'>[FYW]</span>", Motif, fixed = TRUE)) %>% 
    mutate(Motif = gsub("[RHK]", "<span style='color:#9c2175'>[RHK]</span>", Motif, fixed = TRUE))
}

dataset_colors <- c("cTP" = "#9de444", "mTP" = "#e49144", "SP" = "#45e495", "AMP" = "#e44444")

taxonomy_colors <- c("Viridiplantae" = "#ABE188", "Chlorophyta" = "#BFEDB2", "Streptophyta" = "#AAD190",
                     "Unknown" = "grey70", "Metazoa" = "#F7EF99", "Fungi" = "#F1BB87",
                     "Archaea" = "#61A385", "Bacteria" = "#DBA385", "Eukaryota" = "#98D2A4", "Viruses" = "#C66C5D",
                     "Arachnida" = "#ECD0A2", "Insecta" = "#E9B99A", "Lepidosauria" = "#C1D68F", "Mammalia" = "#E6A293",
                     "Arthropoda" = "#EDD4B2", "Chordata" = "#D0A98F")

types <- c("Bigrams without gaps", "Bigrams with gaps", "Trigrams without gaps", 
           "Trigrams with gaps", "Tetragrams without gaps", "Tetragrams with gaps", 
           "Pentagrams without gaps", "Pentagrams with gaps")


###--- Most frequent motifs, full alphabet ---###
df_freq <- read_xlsx(paste0(data_path, "Motifs_results.xlsx"), sheet = "Most frequent") %>% 
  mutate(Motif = sapply(.[["Motif"]], function(i) gsub(".", " _ ", i, fixed = TRUE)))

# Check used cutoff values and set new thresholds if necessary
cutoffs <- df_freq %>% 
  group_by(Type, Dataset, Cutoff) %>% 
  summarise(count = n()) 

# Trigrams with gaps were collected with different cutoffs - set threshold to 0.15 for SP and 0.25 for cTP. Change threshold for AMPs to 0.1
df_freq2 <- df_freq[(which(!(df_freq[["Type"]] == "Trigrams with gaps" & 
                               ((df_freq[["Dataset"]] == "cTP" & df_freq[["Frequency"]] < 0.25) | (df_freq[["Dataset"]] == "SP" & df_freq[["Frequency"]] < 0.2) |
                                  (df_freq[["Dataset"]] == "mTP" & df_freq[["Frequency"]] < 0.125) | (df_freq[["Dataset"]] == "AMP" & df_freq[["Frequency"]] < 0.1))))),]
# Pentagrams with gaps - change cTP threshold to 0.03
df_freq2 <- df_freq2[(which(!(df_freq2[["Type"]] == "Pentagrams with gaps" & df_freq2[["Dataset"]] == "cTP" & df_freq2[["Frequency"]] < 0.03))),]
# Tetragrams with gaps - change cTP and SP thresholds to 0.1
df_freq2 <- df_freq2[(which(!(df_freq2[["Type"]] == "Tetragrams with gaps" & df_freq2[["Dataset"]] %in% c("cTP", "SP") & df_freq2[["Frequency"]] < 0.1))),]
# Tetragrams with gaps - change AMP threshold to 0.05
df_freq2 <- df_freq2[(which(!(df_freq2[["Type"]] == "Tetragrams with gaps" & df_freq2[["Dataset"]] == "AMP" & df_freq2[["Frequency"]] < 0.05))),]

most_frequent_barplots <- lapply(types, function(ith_type) {
  x <- filter(df_freq2, Type == ith_type)
  p <- ggplot(x, aes(x = Frequency, y = reorder_within(Motif, Frequency, Dataset), fill = Dataset)) +
    geom_col() +
    facet_wrap(~Dataset, scales = "free_y", nrow = 1) +
    theme_bw() +
    ggtitle(ith_type) +
    scale_fill_manual("Dataset", values = dataset_colors) +
    scale_y_reordered() +
    theme(legend.position = "none") +
    ylab("Motif")
  ggsave(paste0(data_path, "ngram_results/Most_frequent_barplot_", gsub(" ", "_", ith_type), ".png"),
         p, width = 10, height = 5)
  p
})

most_frequent_barplots_plot <- wrap_plots(most_frequent_barplots, ncol = 2) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D", "E", "F", "G", "H")) +
  plot_layout(design = "
              AAABBBB
              CCCDDDD
              EEEFFFF
              GG#HHHH")
ggsave(paste0(data_path, "ngram_results/Most_frequent_barplot_combined.eps"),
       most_frequent_barplots_plot, width = 20, height = 18)

p1 <- ggplot(df_freq2, aes(x = Type, y = Frequency, color = Dataset)) +
  geom_text(aes(label = Motif), position = position_quasirandom(), size = 3) +
  #geom_quasirandom(size = 0.5) +
  #geom_label_repel(aes(label = Motif), size = 2, label.size = 0.2, label.padding = 0.2) +
  facet_wrap(~Type, scales = "free", nrow = 2) +
  theme_bw() +
  scale_color_manual("Dataset", values = dataset_colors) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(paste0(data_path, "ngram_results/Most_frequent_motifs_all.png"), 
       p1, width = 20, height = 12)

lapply(unique(df_freq2[["Type"]]), function(ith_type) {
  x <- filter(df_freq2, Type == ith_type)
  p <- ggplot(x, aes(x = Type, y = Frequency, color = Dataset)) +
    geom_text(aes(label = Motif), position = position_quasirandom(), size = 3) +
    #geom_quasirandom(size = 0.5) +
    #geom_label_repel(aes(label = Motif), size = 2, label.size = 0.2, label.padding = 0.2) +
    facet_wrap(~Type, scales = "free", nrow = 1) +
    scale_color_manual("Dataset", values = dataset_colors) +
    theme_bw()
  ggsave(paste0(data_path, "ngram_results/Most_frequent_motifs_", gsub(" ", "_", ith_type), ".png"), 
         p, width = 7, height = 10)
  ggsave(paste0(data_path, "ngram_results/Most_frequent_motifs_", gsub(" ", "_", ith_type), ".eps"), 
         p, width = 7, height = 10)
})  



###--- Most frequent motifs, reduced alphabets ---###
df_freq_enc <- read_xlsx(paste0(data_path, "Motifs_results2.xlsx"), sheet = "Most frequent enc") %>% 
  mutate(Motif = sapply(.[["Motif"]], function(i) gsub(".", " _ ", i, fixed = TRUE)))

df_freq_decoded <- lapply(unique(df_freq_enc[["Encoding"]]), function(ith_enc) {
  x <- filter(df_freq_enc, Encoding == ith_enc)
  mutate(x, Motif = sapply(x[["Motif"]], function(i) paste0(decode_motifs(unlist(strsplit(i, "")), ith_enc), collapse = "")))
}) %>% bind_rows()

ggplot(df_freq_decoded, aes(x = Encoding, y = Frequency, fill = Dataset)) +
  geom_violin() +
  facet_wrap(~Type) +
  theme_bw()
ggplot(df_freq_decoded, aes(x = Encoding, y = Frequency, color = Dataset)) +
  geom_quasirandom() +
  facet_wrap(~Type) +
  theme_bw()


# Change pentagrams with gaps thresholds for SP enc7 and enc8 to 0.3
df_freq_decoded2 <- df_freq_decoded[(which(!(df_freq_decoded[["Type"]] == "Pentagrams with gaps" & df_freq_decoded[["Dataset"]] == "SP" & 
                                               df_freq_decoded[["Encoding"]] %in% c("enc7", "enc8") & df_freq_decoded[["Frequency"]] < 0.3))),]
# Change pentagrams with gaps thresholds for cTP enc7 to 0.45
df_freq_decoded2 <- df_freq_decoded2[(which(!(df_freq_decoded2[["Type"]] == "Pentagrams with gaps" & df_freq_decoded2[["Dataset"]] == "cTP" & 
                                                df_freq_decoded2[["Encoding"]] == "enc7" & df_freq_decoded2[["Frequency"]] < 0.45))),]
# Change tetragrams with gaps thresholds for mTP enc8 to 0.3
df_freq_decoded2 <- df_freq_decoded2[(which(!(df_freq_decoded2[["Type"]] == "Tetragrams with gaps" & df_freq_decoded2[["Dataset"]] == "mTP" & 
                                                df_freq_decoded2[["Encoding"]] == "enc8" & df_freq_decoded2[["Frequency"]] < 0.3))),]
# Change tetragrams with gaps thresholds for cTP enc7 to 0.45
df_freq_decoded2 <- df_freq_decoded2[(which(!(df_freq_decoded2[["Type"]] == "Tetragrams with gaps" & df_freq_decoded2[["Dataset"]] == "cTP" & 
                                                df_freq_decoded2[["Encoding"]] == "enc7" & df_freq_decoded2[["Frequency"]] < 0.45))),]
# Change tetragrams with gaps thresholds for SP enc10 to 0.2
df_freq_decoded2 <- df_freq_decoded2[(which(!(df_freq_decoded2[["Type"]] == "Tetragrams with gaps" & df_freq_decoded2[["Dataset"]] == "SP" & 
                                                df_freq_decoded2[["Encoding"]] == "enc10" & df_freq_decoded2[["Frequency"]] < 0.2))),]
# Change trigrams with gaps thresholds for cTP enc7 to 0.75
df_freq_decoded2 <- df_freq_decoded2[(which(!(df_freq_decoded2[["Type"]] == "Trigrams with gaps" & df_freq_decoded2[["Dataset"]] == "cTP" & 
                                                df_freq_decoded2[["Encoding"]] == "enc7" & df_freq_decoded2[["Frequency"]] < 0.75))),]
# Change trigrams with gaps thresholds for cTP enc8 to 0.65
df_freq_decoded2 <- df_freq_decoded2[(which(!(df_freq_decoded2[["Type"]] == "Trigrams with gaps" & df_freq_decoded2[["Dataset"]] == "cTP" & 
                                                df_freq_decoded2[["Encoding"]] == "enc8" & df_freq_decoded2[["Frequency"]] < 0.65))),]
# Change trigrams with gaps thresholds for cTP enc10 to 0.35
df_freq_decoded2 <- df_freq_decoded2[(which(!(df_freq_decoded2[["Type"]] == "Trigrams with gaps" & df_freq_decoded2[["Dataset"]] == "cTP" & 
                                                df_freq_decoded2[["Encoding"]] == "enc10" & df_freq_decoded2[["Frequency"]] < 0.35))),]
# Change trigrams with gaps thresholds for mTP enc7 and enc8 to 0.6
df_freq_decoded2 <- df_freq_decoded2[(which(!(df_freq_decoded2[["Type"]] == "Trigrams with gaps" & df_freq_decoded2[["Dataset"]] == "mTP" & 
                                                df_freq_decoded2[["Encoding"]] %in% c("enc7", "enc8") & df_freq_decoded2[["Frequency"]] < 0.6))),]
# Change trigrams with gaps thresholds for SP enc7 and enc8 to 0.6
df_freq_decoded2 <- df_freq_decoded2[(which(!(df_freq_decoded2[["Type"]] == "Trigrams with gaps" & df_freq_decoded2[["Dataset"]] == "SP" & 
                                                df_freq_decoded2[["Encoding"]] %in% c("enc7", "enc8") & df_freq_decoded2[["Frequency"]] < 0.6))),]
# Change trigrams with gaps thresholds for SP enc10 to 0.5
df_freq_decoded2 <- df_freq_decoded2[(which(!(df_freq_decoded2[["Type"]] == "Trigrams with gaps" & df_freq_decoded2[["Dataset"]] == "SP" & 
                                                df_freq_decoded2[["Encoding"]] == "enc10" & df_freq_decoded2[["Frequency"]] < 0.5))),]



enc_plots <- lapply(types[3:length(types)], function(ith_type) {
  lapply(unique(df_freq_decoded2[["Encoding"]]), function(ith_enc) {
    x <- filter(df_freq_decoded2, Type == ith_type, Encoding == ith_enc) %>% 
      color_encodings()
    p <- ggplot(x, aes(x = Frequency, y = reorder_within(Motif, Frequency, Dataset), fill = Dataset)) +
      geom_col() +
      facet_wrap(~Dataset, scales = "free_y", nrow = 1) +
      theme_bw(base_size = 8) +
      ggtitle(paste0(ith_type, " (encoding: ", ith_enc, ")")) +
      scale_fill_manual("Dataset", values = dataset_colors) +
      ylab("Motif") +
      scale_y_reordered() +
      theme(legend.position = "none",
            axis.text.y = element_markdown())
    #ggsave(paste0(data_path, "ngram_results/Most_frequent_encoded_", gsub(" ", "_", ith_type), "_", ith_enc, ".png"),
    #       p, width = 5*length(unique(x[["Dataset"]])), height = 3+nrow(x)*0.1)
    #ggsave(paste0(data_path, "ngram_results/Most_frequent_encoded_", gsub(" ", "_", ith_type), "_", ith_enc, ".eps"),
    #       p, width = 5*length(unique(x[["Dataset"]])), height = 3+nrow(x)*0.1)
    p
  })
})

enc7_plots <- sapply(enc_plots, function(i) i[1])
enc8_plots <- sapply(enc_plots, function(i) i[2])
enc10_plots <- sapply(enc_plots, function(i) i[3])

enc7_plot <- wrap_plots(enc7_plots, ncol = 2) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D", "E", "F")) +
  plot_layout(design = "
              AAAABBBBBB
              AAAABBBBBB
              AAAABBBBBB
              AAAABBBBBB
              CCCCBBBBBB
              CCCCBBBBBB
              EEEEBBBBBB
              EEEEDDDDDD
              EEEEDDDDDD
              #FF#DDDDDD
              #FF#DDDDDD
              #FF#DDDDDD
              #FF#DDDDDD
              #FF#DDDDDD
              #FF#DDDDDD")
enc7_without_gaps <- wrap_plots(c(enc7_plots[1], enc7_plots[3], enc7_plots[5]), ncol = 1) + 
  plot_annotation(tag_levels = c("A", "B", "C")) +
  plot_layout(design = "
              AAAA
              AAAA
              AAAA
              BBBB
              CCC#
              CCC#")
enc7_with_gaps <- wrap_plots(c(enc7_plots[2], enc7_plots[4], enc7_plots[6]), ncol = 1) + 
  plot_annotation(tag_levels = c("A", "B", "C")) +
  plot_layout(design = "
              AAAAA
              AAAAA
              AAAAA
              AAAAA
              AAAAA
              BBBBB
              BBBBB
              BBBBB
              BBBBB
              BBBBB
              C####
              C####
              C####")
ggsave(paste0(data_path, "ngram_results/Most_frequent_enc7_combined.eps"),
       enc7_plot, width = 20, height = 18)
ggsave(paste0(data_path, "ngram_results/Most_frequent_enc7_nogaps_combined.eps"),
       enc7_without_gaps, width = 8, height = 7)
ggsave(paste0(data_path, "ngram_results/Most_frequent_enc7_gaps_combined.eps"),
       enc7_with_gaps, width = 8, height = 12)

enc8_without_gaps <- wrap_plots(c(enc8_plots[1], enc8_plots[3], enc8_plots[5]), ncol = 1) + 
  plot_annotation(tag_levels = c("A", "B", "C")) +
  plot_layout(design = "
              AAAA
              AAAA
              AAAA
              BBBB
              BBBB
              C###")
enc8_with_gaps <- wrap_plots(c(enc8_plots[2], enc8_plots[4], enc8_plots[6]), ncol = 1) + 
  plot_annotation(tag_levels = c("A", "B", "C")) +
  plot_layout(design = "
              AAAAA
              AAAAA
              AAAAA
              BBBBB
              BBBBB
              BBBBB
              BBBBB
              BBBBB
              CCC##
              CCC##
              CCC##
              CCC##")
ggsave(paste0(data_path, "ngram_results/Most_frequent_enc8_nogaps_combined.eps"),
       enc8_without_gaps, width = 8, height = 7)
ggsave(paste0(data_path, "ngram_results/Most_frequent_enc8_gaps_combined.eps"),
       enc8_with_gaps, width = 8, height = 12)

enc10_without_gaps <- wrap_plots(c(enc10_plots[1], enc10_plots[3], enc10_plots[5]), ncol = 1) + 
  plot_annotation(tag_levels = c("A", "B", "C")) +
  plot_layout(design = "
              AAA
              AAA
              AAA
              AAA
              BB#
              BB#
              BB#
              C##")
enc10_with_gaps <- wrap_plots(c(enc10_plots[2], enc10_plots[4], enc10_plots[6]), ncol = 1) + 
  plot_annotation(tag_levels = c("A", "B", "C")) +
  plot_layout(design = "
              AAAAA
              AAAAA
              AAAAA
              BBBBB
              BBBBB
              BBBBB
              BBBBB
              BBBBB
              BBBBB
              CCC##
              CCC##
              CCC##
              CCC##")
ggsave(paste0(data_path, "ngram_results/Most_frequent_enc10_nogaps_combined.eps"),
       enc10_without_gaps, width = 8, height = 7)
ggsave(paste0(data_path, "ngram_results/Most_frequent_enc10_gaps_combined.eps"),
       enc10_with_gaps, width = 8, height = 12)

###--- Most frequent motifs, taxonomy ---###
df_freq_tax <- read_xlsx(paste0(data_path, "Motifs_results.xlsx"), sheet = "Taxonomy frequent")[,1:10] %>% 
  mutate(Motif = sapply(.[["Motif"]], function(i) gsub(".", " _ ", i, fixed = TRUE))) %>% 
  mutate(`Frequent in` = ifelse(`Frequent in` == "NA", "Unknown", `Frequent in`)) %>% 
  filter(`Frequent in` != "Unknown") %>% 
  mutate(Dataset = case_when(Dataset == "SP" & `Frequent in` %in% c("Viridiplantae", "Metazoa", "Fungi") ~ "SP (kingdom)",
                             Dataset == "SP" & `Frequent in` %in% c("Arthropoda", "Chordata") ~ "SP (phylum)",
                             Dataset == "SP" & `Frequent in` %in% c("Arachnida", "Insecta", "Lepidosauria", "Mammalia") ~ "SP (class)",
                             TRUE ~ Dataset))

# Change bigrams with gaps thresholds for SP Archaea to 0.5
df_freq_tax2 <- df_freq_tax[(which(!(df_freq_tax[["Type"]] == "Bigrams with gaps" & df_freq_tax[["Dataset"]] == "SP" & 
                                       df_freq_tax[["Frequent in"]] == "Archaea" & df_freq_tax[["Freq3"]] < 0.5))),]
# Change pentagrams with gaps thresholds for cTP Chlorophyta to 0.15
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Pentagrams with gaps" & df_freq_tax2[["Dataset"]] == "cTP" & 
                                        df_freq_tax2[["Frequent in"]] == "Chlorophyta" & df_freq_tax2[["Freq2"]] < 0.15))),]
# Change pentagrams with gaps thresholds for cTP Streptophyta to 0.05
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Pentagrams with gaps" & df_freq_tax2[["Dataset"]] == "cTP" & 
                                        df_freq_tax2[["Frequent in"]] == "Streptophyta" & df_freq_tax2[["Freq1"]] < 0.05))),]
# Change pentagrams with gaps thresholds for cTP Viridiplantae to 0.05
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Pentagrams with gaps" & df_freq_tax2[["Dataset"]] == "cTP" & 
                                        df_freq_tax2[["Frequent in"]] == "Viridiplantae" & df_freq_tax2[["Freq1"]] < 0.05))),]
# Change pentagrams with gaps thresholds for mTP Fungi and Metazoa to 0.03
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Pentagrams with gaps" & df_freq_tax2[["Dataset"]] == "mTP" & 
                                        df_freq_tax2[["Frequent in"]] == "Fungi" & df_freq_tax2[["Freq3"]] < 0.03))),]
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Pentagrams with gaps" & df_freq_tax2[["Dataset"]] == "mTP" & 
                                        df_freq_tax2[["Frequent in"]] == "Metazoa" & df_freq_tax2[["Freq2"]] < 0.03))),]
# Change pentagrams with gaps thresholds for mTP Viridiplantae to 0.05
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Pentagrams with gaps" & df_freq_tax2[["Dataset"]] == "mTP" & 
                                        df_freq_tax2[["Frequent in"]] == "Viridiplantae" & df_freq_tax2[["Freq1"]] < 0.05))),]
# Change pentagrams with gaps thresholds for SP Archaea to 0.2
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Pentagrams with gaps" & df_freq_tax2[["Dataset"]] == "SP" & 
                                        df_freq_tax2[["Frequent in"]] == "Archaea" & df_freq_tax2[["Freq3"]] < 0.2))),]
# Change pentagrams with gaps thresholds for SP Eukaryota to 0.05
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Pentagrams with gaps" & df_freq_tax2[["Dataset"]] == "SP" & 
                                        df_freq_tax2[["Frequent in"]] == "Eukaryota" & df_freq_tax2[["Freq1"]] < 0.05))),]
# Change pentagrams with gaps thresholds for SP Viruses to 0.1
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Pentagrams with gaps" & df_freq_tax2[["Dataset"]] == "SP" & 
                                        df_freq_tax2[["Frequent in"]] == "Viruses" & df_freq_tax2[["Freq4"]] < 0.1))),]
# Change pentagrams with gaps thresholds for SP Archaea to 0.2
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Pentagrams with gaps" & df_freq_tax2[["Dataset"]] == "SP" & 
                                        df_freq_tax2[["Frequent in"]] == "Archaea" & df_freq_tax2[["Freq3"]] < 0.2))),]
# Change tetragrams with gaps thresholds for cTP Chlorophyta to 0.25
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Tetragrams with gaps" & df_freq_tax2[["Dataset"]] == "cTP" & 
                                        df_freq_tax2[["Frequent in"]] == "Chlorophyta" & df_freq_tax2[["Freq2"]] < 0.25))),]
# Change tetragrams with gaps thresholds for cTP Streptophyta and Viridiplantae to 0.1
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Tetragrams with gaps" & df_freq_tax2[["Dataset"]] == "cTP" & 
                                        df_freq_tax2[["Frequent in"]] %in% c("Viridiplantae", "Streptophyta") & df_freq_tax2[["Freq1"]] < 0.1))),]
# Change tetragrams with gaps thresholds for mTP Metazoa to 0.075
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Tetragrams with gaps" & df_freq_tax2[["Dataset"]] == "mTP" & 
                                        df_freq_tax2[["Frequent in"]] == "Metazoa" & df_freq_tax2[["Freq2"]] < 0.075))),]
# Change tetragrams with gaps thresholds for mTP Viridiplantae to 0.075
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Tetragrams with gaps" & df_freq_tax2[["Dataset"]] == "mTP" & 
                                        df_freq_tax2[["Frequent in"]] == "Viridiplantae" & df_freq_tax2[["Freq1"]] < 0.075))),]
# Change tetragrams with gaps thresholds for SP Archaea to 0.25
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Tetragrams with gaps" & df_freq_tax2[["Dataset"]] == "SP" & 
                                        df_freq_tax2[["Frequent in"]] == "Archaea" & df_freq_tax2[["Freq3"]] < 0.25))),]
# Change tetragrams with gaps thresholds for SP Eukaryota to 0.1
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Tetragrams with gaps" & df_freq_tax2[["Dataset"]] == "SP" & 
                                        df_freq_tax2[["Frequent in"]] == "Eukaryota" & df_freq_tax2[["Freq1"]] < 0.1))),]
# Change tetragrams with gaps thresholds for SP Viruses to 0.15
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Tetragrams with gaps" & df_freq_tax2[["Dataset"]] == "SP" & 
                                        df_freq_tax2[["Frequent in"]] == "Viruses" & df_freq_tax2[["Freq4"]] < 0.15))),]
# Change trigrams with gaps thresholds for cTP Chlorophyta to 0.4
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Trigrams with gaps" & df_freq_tax2[["Dataset"]] == "cTP" & 
                                        df_freq_tax2[["Frequent in"]] == "Chlorophyta" & df_freq_tax2[["Freq2"]] < 0.4))),]
# Change trigrams with gaps thresholds for cTP Streptophyta to 0.3
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Trigrams with gaps" & df_freq_tax2[["Dataset"]] == "cTP" & 
                                        df_freq_tax2[["Frequent in"]] %in% c("Viridiplantae", "Streptophyta") & df_freq_tax2[["Freq1"]] < 0.3))),]
# Change trigrams with gaps thresholds for mTP Metazoa to 0.15
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Trigrams with gaps" & df_freq_tax2[["Dataset"]] == "mTP" & 
                                        df_freq_tax2[["Frequent in"]] == "Metazoa" & df_freq_tax2[["Freq2"]] < 0.15))),]
# Change trigrams with gaps thresholds for mTP Viridiplantae to 0.2
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Trigrams with gaps" & df_freq_tax2[["Dataset"]] == "mTP" & 
                                        df_freq_tax2[["Frequent in"]] == "Viridiplantae" & df_freq_tax2[["Freq1"]] < 0.2))),]
# Change trigrams with gaps thresholds for SP Bacteria to 0.2
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Trigrams with gaps" & df_freq_tax2[["Dataset"]] == "SP" & 
                                        df_freq_tax2[["Frequent in"]] == "Bacteria" & df_freq_tax2[["Freq2"]] < 0.2))),]
# Change trigrams with gaps thresholds for SP Eukaryota to 0.3
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Trigrams with gaps" & df_freq_tax2[["Dataset"]] == "SP" & 
                                        df_freq_tax2[["Frequent in"]] == "Eukaryota" & df_freq_tax2[["Freq1"]] < 0.3))),]
# Change trigrams with gaps thresholds for SP Archaea to 0.35
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Trigrams with gaps" & df_freq_tax2[["Dataset"]] == "SP" & 
                                        df_freq_tax2[["Frequent in"]] == "Archaea" & df_freq_tax2[["Freq3"]] < 0.35))),]
# Change trigrams with gaps thresholds for SP Viruses to 0.25
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Trigrams with gaps" & df_freq_tax2[["Dataset"]] == "SP" & 
                                        df_freq_tax2[["Frequent in"]] == "Viruses" & df_freq_tax2[["Freq4"]] < 0.25))),]
# Change pentagrams with gaps thresholds for SP Arachnida to 0.07
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Pentagrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (class)" & 
                                        df_freq_tax2[["Frequent in"]] == "Arachnida" & df_freq_tax2[["Freq1"]] < 0.07))),]
# Change pentagrams with gaps thresholds for SP Fungi to 0.03
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Pentagrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (kingdom)" & 
                                        df_freq_tax2[["Frequent in"]] == "Fungi" & df_freq_tax2[["Freq3"]] < 0.03))),]
# Change tetragrams with gaps thresholds for SP Arachnida to 0.1
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Tetragrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (class)" & 
                                        df_freq_tax2[["Frequent in"]] == "Arachnida" & df_freq_tax2[["Freq1"]] < 0.1))),]
# Change tetragrams with gaps thresholds for SP Insecta to 0.055
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Tetragrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (class)" & 
                                        df_freq_tax2[["Frequent in"]] == "Insecta" & df_freq_tax2[["Freq2"]] < 0.055))),]
# Change tetragrams with gaps thresholds for SP Lepidosauria to 0.2
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Tetragrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (class)" & 
                                        df_freq_tax2[["Frequent in"]] == "Lepidosauria" & df_freq_tax2[["Freq3"]] < 0.2))),]
# Change tetragrams with gaps thresholds for SP Mammalia to 0.2
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Tetragrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (class)" & 
                                        df_freq_tax2[["Frequent in"]] == "Mammalia" & df_freq_tax2[["Freq4"]] < 0.2))),]
# Change tetragrams with gaps thresholds for SP Fungi to 0.06
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Tetragrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (kingdom)" & 
                                        df_freq_tax2[["Frequent in"]] == "Fungi" & df_freq_tax2[["Freq3"]] < 0.06))),]
# Change tetragrams with gaps thresholds for SP Metazoa to 0.125
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Tetragrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (kingdom)" & 
                                        df_freq_tax2[["Frequent in"]] == "Metazoa" & df_freq_tax2[["Freq2"]] < 0.125))),]
# Change tetragrams with gaps thresholds for SP Viridiplantae to 0.0625
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Tetragrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (kingdom)" & 
                                        df_freq_tax2[["Frequent in"]] == "Viridiplantae" & df_freq_tax2[["Freq1"]] < 0.0625))),]
# Change tetragrams with gaps thresholds for SP Chordata to 0.1
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Tetragrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (phylum)" & 
                                        df_freq_tax2[["Frequent in"]] == "Chordata" & df_freq_tax2[["Freq2"]] < 0.1))),]
# Change trigrams with gaps thresholds for SP Arachnida to 0.2
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Trigrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (class)" & 
                                        df_freq_tax2[["Frequent in"]] == "Arachnida" & df_freq_tax2[["Freq1"]] < 0.2))),]
# Change trigrams with gaps thresholds for SP Lepidosauria to 0.3
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Trigrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (class)" & 
                                        df_freq_tax2[["Frequent in"]] == "Lepidosauria" & df_freq_tax2[["Freq3"]] < 0.3))),]
# Change tetragrams with gaps thresholds for SP Mammalia to 0.25
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Trigrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (class)" & 
                                        df_freq_tax2[["Frequent in"]] == "Mammalia" & df_freq_tax2[["Freq4"]] < 0.25))),]
# Change trigrams with gaps thresholds for SP Fungi to 0.2
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Trigrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (kingdom)" & 
                                        df_freq_tax2[["Frequent in"]] == "Fungi" & df_freq_tax2[["Freq3"]] < 0.2))),]
# Change trigrams with gaps thresholds for SP Metazoa to 0.25
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Trigrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (kingdom)" & 
                                        df_freq_tax2[["Frequent in"]] == "Metazoa" & df_freq_tax2[["Freq2"]] < 0.25))),]
# Change trigrams with gaps thresholds for SP Viridiplantae to 0.175
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Trigrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (kingdom)" & 
                                        df_freq_tax2[["Frequent in"]] == "Viridiplantae" & df_freq_tax2[["Freq1"]] < 0.175))),]
# Change trigrams with gaps thresholds for SP Chordata to 0.2
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Type"]] == "Trigrams with gaps" & df_freq_tax2[["Dataset"]] == "SP (phylum)" & 
                                        df_freq_tax2[["Frequent in"]] == "Chordata" & df_freq_tax2[["Freq2"]] < 0.2))),]

# Remove Vidiriplantae/NA from cTP
df_freq_tax2 <- df_freq_tax2[(which(!(df_freq_tax2[["Dataset"]] == "cTP" & df_freq_tax2[["Frequent in"]] == "Viridiplantae"))),]

saveRDS(df_freq_tax2, "./data/df_freq_tax2.rds", compress = "xz")

tax_plots <- lapply(unique(df_freq_tax2[["Dataset"]]), function(ith_set) {
  lapply(types, function(ith_type) {
    x <- filter(df_freq_tax2, Dataset == ith_set, Type == ith_type) %>% 
      select(-Diff) %>% 
      pivot_longer(c("Freq1", "Freq2", "Freq3", "Freq4"), names_to = "Group", values_to = "Frequency") %>% 
      filter(!is.na(Frequency))
    
    if(ith_set == "cTP") {
      y <- mutate(x, Group = case_when(Group == "Freq1" & `Frequent in` %in% c("Viridiplantae", "Unknown") ~ "Viridiplantae",
                                       Group == "Freq2" & `Frequent in` %in% c("Viridiplantae", "Unknown") ~ "Unknown",
                                       Group == "Freq1" & `Frequent in` %in% c("Streptophyta", "Chlorophyta") ~ "Streptophyta",
                                       Group == "Freq2" & `Frequent in` %in% c("Streptophyta", "Chlorophyta") ~ "Chlorophyta")) 
    } else if(ith_set == "mTP") {
      y <- mutate(x, Group = case_when(Group == "Freq1" ~ "Viridiplantae",
                                       Group == "Freq2" ~ "Metazoa",
                                       Group == "Freq3" ~ "Fungi",
                                       Group == "Freq4" ~ "Unknown")) %>% 
        filter(Group != "Unknown")
    } else if(ith_set == "SP") {
      y <- mutate(x, Group = case_when(Group == "Freq1" ~ "Eukaryota",
                                       Group == "Freq2" ~ "Bacteria",
                                       Group == "Freq3" ~ "Archaea",
                                       Group == "Freq4" ~ "Viruses")) 
    } else if(ith_set == "SP (kingdom)") {
      y <- mutate(x, Group = case_when(Group == "Freq1" ~ "Viridiplantae",
                                       Group == "Freq2" ~ "Metazoa",
                                       Group == "Freq3" ~ "Fungi")) 
    } else if(ith_set == "SP (phylum)") {
      y <- mutate(x, Group = case_when(Group == "Freq1" ~ "Arthropoda",
                                       Group == "Freq2" ~ "Chordata")) 
    } else if(ith_set == "SP (class)") {
      y <- mutate(x, Group = case_when(Group == "Freq1" ~ "Arachnida",
                                       Group == "Freq2" ~ "Insecta",
                                       Group == "Freq3" ~ "Lepidosauria",
                                       Group == "Freq4" ~ "Mammalia")) 
    }
    p <- ggplot(y, aes(y = Motif, x = Frequency, fill = Group)) +
      geom_col(position = position_dodge()) +
      facet_wrap(~`Frequent in`, scales = "free", nrow = 1) +
      theme_bw(base_size = 6) +
      ggtitle(paste0(ith_type, ", ", ith_set)) +
      scale_fill_manual("Taxonomic group", values = taxonomy_colors[names(taxonomy_colors) %in% unique(y[["Group"]])])
    ggsave(paste0(data_path, "ngram_results/Taxonomy_all_", gsub(" ", "_", ith_type), "_", gsub(" ", "_", ith_set), ".png"),
           p, width = 10, height = 2+nrow(y)*0.02, limitsize = FALSE)
    ggsave(paste0(data_path, "ngram_results/Taxonomy_all_", gsub(" ", "_", ith_type), "_", gsub(" ", "_", ith_set), ".eps"),
           p, width = 8, height = 2+nrow(y)*0.02, limitsize = FALSE)
    p
  })
})

tax_ctp <- wrap_plots(tax_plots[[1]], ncol = 2) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D", "E", "F", "G", "H")) +
  plot_layout(design = "
              AABB
              AABB
              AABB
              CCBB
              CCBB
              DDEE
              DDEE
              FFGG
              FFGG
              FFHH
              FFHH", guides = 'collect') & theme(legend.position = "bottom")
ggsave(paste0(data_path, "ngram_results/Taxonomy_all_cTP_combined.eps"),
       tax_ctp, width = 8, height = 12)

tax_mtp <- wrap_plots(tax_plots[[2]], ncol = 2) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D", "E", "F", "G", "H")) +
  plot_layout(design = "
              AABB
              AABB
              CCBB
              DDBB
              DDEE
              DDEE
              DDEE
              FFGG
              FFGG
              FFHH
              FFHH", guides = 'collect') & theme(legend.position = "bottom")
ggsave(paste0(data_path, "ngram_results/Taxonomy_all_mTP_combined.eps"),
       tax_mtp, width = 9, height = 12)

tax_sp <- wrap_plots(tax_plots[[3]], ncol = 2) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D", "E", "F", "G", "H")) +
  plot_layout(design = "
              AABBB
              AABBB
              AABBB
              CCBBB
              CCDDD
              EEDDD
              EEFFF
              GGFFF
              GGFFF
              GGFFF
              GGFFF
              HHFFF
              HHFFF
              HHFFF", guides = 'collect') & theme(legend.position = "bottom")
ggsave(paste0(data_path, "ngram_results/Taxonomy_all_SP_superkingdom_combined.eps"),
       tax_sp, width = 12, height = 12)

tax_sp_k <- wrap_plots(tax_plots[[4]], ncol = 2) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D", "E", "F", "G", "H")) +
  plot_layout(design = "
              AABB
              CCBB
              CCBB
              DDEE
              DDEE
              DDFF
              GGFF
              HHFF
              HHFF
              HHFF", guides = 'collect') & theme(legend.position = "bottom")
ggsave(paste0(data_path, "ngram_results/Taxonomy_all_SP_kingdom_combined.eps"),
       tax_sp_k, width = 10, height = 11)

tax_sp_p <- wrap_plots(tax_plots[[5]], ncol = 2) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D", "E", "F", "G", "H")) +
  plot_layout(design = "
              AABB
              CCBB
              CCDD
              EEDD
              EEDD
              EEDD
              FFG#
              FFG#
              FFG#
              FFHH
              FFHH
              FFHH
              FFHH", guides = 'collect') & theme(legend.position = "bottom")
ggsave(paste0(data_path, "ngram_results/Taxonomy_all_SP_phylum_combined.eps"),
       tax_sp_p, width = 9, height = 12)

tax_sp_c <- wrap_plots(tax_plots[[6]], ncol = 2) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D", "E", "F", "G", "H")) +
  plot_layout(design = "
              AABB
              AABB
              CCDD
              CCDD
              EEDD
              EEDD
              EEGG
              FFHH
              FFHH
              FFHH
              FFHH", guides = 'collect') & theme(legend.position = "bottom")
ggsave(paste0(data_path, "ngram_results/Taxonomy_all_SP_class_combined.eps"),
       tax_sp_c, width = 12, height = 13)


tax_group_plots <- lapply(unique(df_freq_tax2[["Dataset"]]), function(ith_set) {
  lapply(types, function(ith_type) {
    x <- filter(df_freq_tax2, Dataset == ith_set, Type == ith_type) %>% 
      select(-Diff) %>% 
      pivot_longer(c("Freq1", "Freq2", "Freq3", "Freq4"), names_to = "Group", values_to = "Frequency") %>% 
      filter(!is.na(Frequency))
    
    if(ith_set == "cTP") {
      y <- mutate(x, Group = case_when(Group == "Freq1" & `Frequent in` %in% c("Viridiplantae", "Unknown") ~ "Viridiplantae",
                                       Group == "Freq2" & `Frequent in` %in% c("Viridiplantae", "Unknown") ~ "Unknown",
                                       Group == "Freq1" & `Frequent in` %in% c("Streptophyta", "Chlorophyta") ~ "Streptophyta",
                                       Group == "Freq2" & `Frequent in` %in% c("Streptophyta", "Chlorophyta") ~ "Chlorophyta")) 
    } else if(ith_set == "mTP") {
      y <- mutate(x, Group = case_when(Group == "Freq1" ~ "Viridiplantae",
                                       Group == "Freq2" ~ "Metazoa",
                                       Group == "Freq3" ~ "Fungi",
                                       Group == "Freq4" ~ "Unknown")) %>% 
        filter(Group != "Unknown")
    } else if(ith_set == "SP") {
      y <- mutate(x, Group = case_when(Group == "Freq1" ~ "Eukaryota",
                                       Group == "Freq2" ~ "Bacteria",
                                       Group == "Freq3" ~ "Archaea",
                                       Group == "Freq4" ~ "Viruses")) 
    } else if(ith_set == "SP (kingdom)") {
      y <- mutate(x, Group = case_when(Group == "Freq1" ~ "Viridiplantae",
                                       Group == "Freq2" ~ "Metazoa",
                                       Group == "Freq3" ~ "Fungi")) 
    } else if(ith_set == "SP (phylum)") {
      y <- mutate(x, Group = case_when(Group == "Freq1" ~ "Arthropoda",
                                       Group == "Freq2" ~ "Chordata")) 
    } else if(ith_set == "SP (class)") {
      y <- mutate(x, Group = case_when(Group == "Freq1" ~ "Arachnida",
                                       Group == "Freq2" ~ "Insecta",
                                       Group == "Freq3" ~ "Lepidosauria",
                                       Group == "Freq4" ~ "Mammalia")) 
    }
    p <- filter(y, Group == `Frequent in`) %>% 
      ggplot(aes(y = reorder_within(Motif, Frequency, Group), x = Frequency, fill = Group)) +
      geom_col(position = position_dodge()) +
      scale_y_reordered() +
      ylab("Motif") +
      facet_wrap(~`Frequent in`, scales = "free", nrow = 1) +
      theme_bw(base_size = 6) +
      ggtitle(paste0(ith_type, ", ", ith_set)) +
      scale_fill_manual("Taxonomic group", values = taxonomy_colors[names(taxonomy_colors) %in% unique(y[["Group"]])])
    ggsave(paste0(data_path, "ngram_results/Taxonomy_by_group_", gsub(" ", "_", ith_type), "_", gsub(" ", "_", ith_set), ".png"),
           p, width = 10, height = 2+nrow(y)*0.02, limitsize = FALSE)
    ggsave(paste0(data_path, "ngram_results/Taxonomy_by_group_", gsub(" ", "_", ith_type), "_", gsub(" ", "_", ith_set), ".eps"),
           p, width = 10, height = 2+nrow(y)*0.02, limitsize = FALSE)
    p
  })
})

tax_groups_ctp <- wrap_plots(tax_group_plots[[1]], ncol = 2) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D", "E", "F", "G", "H")) +
  plot_layout(design = "
              AABB
              AABB
              AABB
              AABB
              AABB
              CCBB
              CCBB
              CCDD
              EEDD
              EEDD
              EEDD
              EEDD
              FFGG
              FFGG
              FFGG
              FFGG
              FFHH
              FFHH
              FFHH", guides = 'collect') & theme(legend.position = "none")
ggsave(paste0(data_path, "ngram_results/Taxonomy_groups_cTP_combined.eps"),
       tax_groups_ctp, width = 8, height = 12)

tax_groups_mtp <- wrap_plots(tax_group_plots[[2]], ncol = 2) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D", "E", "F", "G", "H")) +
  plot_layout(design = "
              AABB
              AABB
              CCBB
              DDBB
              DDEE
              DDEE
              DDEE
              FFGG
              FFGG
              FFHH
              FFHH", guides = 'collect') & theme(legend.position = "none")
ggsave(paste0(data_path, "ngram_results/Taxonomy_groups_mTP_combined.eps"),
       tax_groups_mtp, width = 10, height = 12)

tax_groups_sp <- wrap_plots(tax_group_plots[[3]], ncol = 2) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D", "E", "F", "G", "H")) +
  plot_layout(design = "
              AABB
              CCBB
              DDEE
              FFGG
              FFGG
              FFGG
              FFHH
              FFHH", guides = 'collect') & theme(legend.position = "none")
ggsave(paste0(data_path, "ngram_results/Taxonomy_groups_SP_superkingdom_combined.eps"),
       tax_groups_sp, width = 12, height = 12)

tax_groups_sp_k <- wrap_plots(tax_group_plots[[4]], ncol = 2) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D", "E", "F", "G", "H")) +
  plot_layout(design = "
              AABB
              CCBB
              DDEE
              DDEE
              FFGG
              FFHH
              FFHH
              FFHH", guides = 'collect') & theme(legend.position = "none")
ggsave(paste0(data_path, "ngram_results/Taxonomy_groups_SP_kingdom_combined.eps"),
       tax_groups_sp_k, width = 9, height = 11)

tax_groups_sp_p <- wrap_plots(tax_group_plots[[5]], ncol = 2) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D", "E", "F", "G", "H")) +
  plot_layout(design = "
              AABB
              CCBB
              DDEE
              DDEE
              DDEE
              DDFF
              G#FF
              G#FF
              HHFF
              HHFF", guides = 'collect') & theme(legend.position = "none")
ggsave(paste0(data_path, "ngram_results/Taxonomy_groups_SP_phylum_combined.eps"),
       tax_groups_sp_p, width = 9, height = 14)

tax_groups_sp_c <- wrap_plots(tax_group_plots[[6]], ncol = 2) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D", "E", "F", "G", "H")) +
  plot_layout(design = "
              AABB
              CCBB
              CCDD
              EEDD
              EEDD
              FFGG
              FFHH
              FFHH
              FFHH
              FFHH
              ", guides = 'collect') & theme(legend.position = "none")
ggsave(paste0(data_path, "ngram_results/Taxonomy_groups_SP_class_combined.eps"),
       tax_groups_sp_c, width = 11, height = 13)


###--- Common motifs, taxonomy ---###
df_common_tax <- read_xlsx(paste0(data_path, "Motifs_results.xlsx"), sheet = "Taxonomy common") %>% 
  mutate(Motif = sapply(.[["Motif"]], function(i) gsub(".", " _ ", i, fixed = TRUE))) %>% 
  # filter cTP to include only motifs with frequency >= 0.1 in Streptophyta & Chlorophyta or frequency >= 0.2 in Viridiplantae & NA
  # filter SP to include only motifs with frequency >= 0.1 in Eukaryota & Archaea
  filter((Groups == "Viridiplantae, NA" & Freq1 >= 0.2 & Freq2 >= 0.2) |
           (Groups == "Streptophyta, Chlorophyta" & Freq1 >= 0.15 & Freq2 >= 0.15) |
           (Groups == "Eukaryota, Archaea" & Freq1 >= 0.1 & Freq3 >= 0.1) |
           Dataset == "mTP" |
           Dataset == "SP" & Groups %in% c("Eukaryota, Bacteria, Archaea, Viruses", "Eukaryota, Bacteria, Archaea",
                                           "Eukaryota, Bacteria", "Eukaryota, Bacteria, Viruses", 
                                           "Eukaryota, Viruses", "Eukaryota, Archaea, Viruses")) 

lapply(unique(df_common_tax[["Dataset"]]), function(ith_set) {
  x <- filter(df_common_tax, Dataset == ith_set) %>% 
    pivot_longer(c("Freq1", "Freq2", "Freq3", "Freq4"), names_to = "Group", values_to = "Frequency") %>% 
    filter(!is.na(Frequency))
  
  if(ith_set == "cTP") {
    y <- mutate(x, Group = case_when(Group == "Freq1" & Groups == "Viridiplantae, NA" ~ "Viridiplantae",
                                     Group == "Freq2" & Groups == "Viridiplantae, NA" ~ "Unknown",
                                     Group == "Freq1" & Groups == "Streptophyta, Chlorophyta" ~ "Streptophyta",
                                     Group == "Freq2" & Groups == "Streptophyta, Chlorophyta" ~ "Chlorophyta"))
  } else if(ith_set == "mTP") {
    y <- mutate(x, Group = case_when(Group == "Freq1" ~ "Viridiplantae",
                                     Group == "Freq2" ~ "Metazoa",
                                     Group == "Freq3" ~ "Fungi",
                                     Group == "Freq4" ~ "Unknown"))  %>% 
      filter(Group != "Unknown")
  } else if(ith_set == "SP") {
    y <- mutate(x, Group = case_when(Group == "Freq1" ~ "Eukaryota",
                                     Group == "Freq2" ~ "Bacteria",
                                     Group == "Freq3" ~ "Archaea",
                                     Group == "Freq4" ~ "Viruses")) 
  } 
  p <- ggplot(y, aes(x = Frequency, y = Motif, fill = Group)) +
    geom_col(position = position_dodge()) +
    facet_wrap(~Groups, scales = "free", nrow = 1) +
    theme_bw() +
    ggtitle(ith_set) +
    scale_fill_manual("Taxonomic group", values = taxonomy_colors[names(taxonomy_colors) %in% unique(y[["Group"]])])
  ggsave(paste0(data_path, "ngram_results/Taxonomy_common_", ith_set, ".png"),
         p, width = 4*length(unique(y[["Groups"]])), height = 4+nrow(y)*0.02, limitsize = FALSE)
  ggsave(paste0(data_path, "ngram_results/Taxonomy_common_", ith_set, ".eps"),
         p, width = 4*length(unique(y[["Groups"]])), height = 4+nrow(y)*0.02, limitsize = FALSE)
})




###--- Differentiating motifs ---###
df_diff <- read_xlsx(paste0(data_path, "Motifs_results.xlsx"), sheet = "Differences") %>% 
  mutate(Motif = sapply(.[["Motif"]], function(i) gsub(".", " _ ", i, fixed = TRUE))) %>% 
  filter(!(Type == "Tetragrams with gaps" & FreqDiff < 0.075 & Comparison == "AMP/cTP")) %>% 
  filter(!(Type == "Trigrams with gaps" & FreqDiff < 0.15))

lapply(unique(df_diff[["Type"]]), function(ith_type) {
  x <- filter(df_diff, Type == ith_type) %>% 
    select(-c("quipt_pval", "Cutoff")) %>% 
    pivot_longer(c(Freq1, Freq2), names_to = "Set", values_to = "Frequency")
  p <- ggplot(x, aes(y = reorder_within(Motif, FreqDiff, Type), x = Frequency, fill = Set)) +
    geom_col(position = position_dodge()) +
    facet_wrap(~Comparison, scales = "free", nrow = 1) +
    theme_bw() +
    scale_y_reordered() +
    ggtitle(ith_type) +
    ylab("Motif")
  ggsave(paste0(data_path, "ngram_results/Differences_", ith_type, ".png"),
         p, width = 4*length(unique(x[["Comparison"]])), height = 4+nrow(x)*0.02, limitsize = FALSE)
})



###--- Taxonomic groups distribution ---###
library(targets)
tax_dat <- read_xlsx(paste0(data_path, "Presequence_uniprot_entries_proper.xlsx"))
tar_load(datasets_list)

tax_plot_dat <- tax_dat %>% 
  mutate(Dataset = case_when(Entry %in% names(datasets_list[["cTP experimentally verified presequence"]]) ~ "cTP",
                             Entry %in% names(datasets_list[["mTP experimentally verified presequence"]]) ~ "mTP",
                             Entry %in% names(datasets_list[["SP experimentally verified presequence"]]) ~ "SP")) %>% 
  filter(!is.na(Dataset))

tax_seq_numbers <- tax_plot_dat %>% 
  group_by(Dataset, superkingdom, kingdom, phylum, class) %>% 
  summarise(count = n()) 
write.csv(tax_seq_numbers, "./data/tax_seq_numbers.csv", row.names = FALSE)

tax_org_numbers <- tax_plot_dat %>% 
  group_by(Dataset, superkingdom, kingdom, phylum, class) %>% 
  summarise(count = length(unique(Organism)))
write.csv(tax_org_numbers , "./data/tax_org_numbers.csv", row.names = FALSE)

write.csv(tax_plot_dat, "./data/tax_plot_dat.csv", row.names = FALSE)

p_ctp <- tax_plot_dat %>% 
  filter(Dataset == "cTP") %>% 
  group_by(Dataset, phylum) %>%
  summarise(`Number of sequences` = n()) %>% 
  ggplot(aes(y = reorder(phylum, `Number of sequences`, decreasing = TRUE), x = `Number of sequences`, fill = Dataset)) +
  geom_col() +
  theme_bw(base_size = 8) +
  geom_text(aes(label = `Number of sequences`), hjust = -0.25, size = 2) +
  scale_fill_manual(values = c("cTP" = "#9de444")) +
  theme(legend.position = "none") +
  ylab("Phylum") +
  ggtitle("Number of cTP sequences in each phylum")
ggsave(paste0(data_path, "ngram_results/taxonomy_numbers_of_sequences/cTP_phylum.png"), p_ctp, width = 6, height = 4)

p_mtp <- tax_plot_dat %>% 
  filter(Dataset == "mTP") %>% 
  group_by(Dataset, kingdom) %>%
  summarise(`Number of sequences` = n()) %>% 
  ggplot(aes(y = reorder(kingdom, `Number of sequences`, decreasing = TRUE), x = `Number of sequences`, fill = Dataset)) +
  geom_col() +
  theme_bw(base_size = 8) +
  geom_text(aes(label = `Number of sequences`), hjust = -0.25, size = 2) +
  scale_fill_manual(values = c("mTP" = "#e49144")) +
  theme(legend.position = "none") +
  ylab("Kingdom") +
  ggtitle("Number of mTP sequences in each kingdom")
ggsave(paste0(data_path, "ngram_results/taxonomy_numbers_of_sequences/mTP_kingdom.png"), p_mtp, width = 6, height = 4)

p_sp <- tax_plot_dat %>% 
  filter(Dataset == "SP") %>% 
  group_by(Dataset, superkingdom) %>%
  summarise(`Number of sequences` = n()) %>% 
  ggplot(aes(y = reorder(superkingdom, `Number of sequences`, decreasing = TRUE), x = `Number of sequences`, fill = Dataset)) +
  geom_col() +
  theme_bw(base_size = 8) +
  geom_text(aes(label = `Number of sequences`), hjust = -0.25, size = 2) +
  scale_fill_manual(values = c("SP" = "#45e495")) +
  theme(legend.position = "none") +
  ylab("Superkingdom") +
  ggtitle("Number of SP sequences in each superkingdom")
ggsave(paste0(data_path, "ngram_results/taxonomy_numbers_of_sequences/SP_superkingdom.png"), p_sp, width = 6, height = 4)


