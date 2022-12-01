library(dplyr)
library(ggplot2)
library(readxl)
library(ggrepel)
library(ggbeeswarm)
library(biogram)
library(tidyr)

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


###--- Most frequent motifs, full alphabet ---###
df_freq <- read_xlsx(paste0(data_path, "Motifs_results.xlsx"), sheet = "Most frequent") %>% 
  mutate(Motif = sapply(.[["Motif"]], function(i) gsub(".", " _ ", i, fixed = TRUE)))

# Check used cutoff values and set new thresholds if necessary
cutoffs <- df_freq %>% 
  group_by(Type, Dataset, Cutoff) %>% 
  summarise(count = n()) 

# Trigrams with gaps were collected with different cutoffs - set threshold to 0.15 for SP and 0.2 for cTP
df_freq2 <- df_freq[(which(!(df_freq[["Type"]] == "Trigrams with gaps" & 
                               ((df_freq[["Dataset"]] == "cTP" & df_freq[["Frequency"]] < 0.2) | (df_freq[["Dataset"]] == "SP" & df_freq[["Frequency"]] < 0.15))))),]

# x <- filter(df_freq2, Type == "Trigrams with gaps")

# ggplot(x, aes(x = Dataset, y = Frequency)) +
#   geom_quasirandom() +
#   facet_wrap(~Dataset, scales = "free", nrow = 1) +
#   theme_bw() +
#   geom_label_repel(aes(label = Motif))

# ggplot(x, aes(y = reorder(Motif, Frequency), x = Frequency)) +
#   geom_col() +
#   facet_wrap(~Dataset, scales = "free", nrow = 1) +
#   theme_bw() 

lapply(unique(df_freq2[["Type"]]), function(ith_type) {
  x <- filter(df_freq2, Type == ith_type)
  p <- lapply(unique(x[["Dataset"]]), function(i) {
    filter(x, Dataset == i) %>% 
      arrange(desc(Frequency)) %>% 
      head(n = 15)
  }) %>% bind_rows() %>% 
    ggplot(aes(x = Frequency, y = reorder(Motif, Frequency))) +
    geom_col() +
    facet_wrap(~Dataset, scales = "free_y", nrow = 1) +
    theme_bw() +
    ggtitle(ith_type)
  ggsave(paste0(data_path, "ngram_results/Most_frequent_15best_", gsub(" ", "_", ith_type), ".png"), 
         p, width = 10, height = 5)
})

p1 <- ggplot(df_freq2, aes(x = Type, y = Frequency, color = Dataset)) +
  geom_text(aes(label = Motif), position = position_quasirandom(), size = 3) +
  #geom_quasirandom(size = 0.5) +
  #geom_label_repel(aes(label = Motif), size = 2, label.size = 0.2, label.padding = 0.2) +
  facet_wrap(~Type, scales = "free", nrow = 1) +
  theme_bw()
ggsave(paste0(data_path, "ngram_results/Most_frequent_motifs.png"), 
       p1, width = 20, height = 12)

lapply(unique(df_freq2[["Type"]]), function(ith_type) {
  x <- filter(df_freq2, Type == ith_type)
  p <- lapply(unique(x[["Dataset"]]), function(i) {
    filter(x, Dataset == i) %>% 
      arrange(desc(Frequency)) %>% 
      head(n = 15)
  }) %>% bind_rows() %>% 
    ggplot(aes(x = Type, y = Frequency, color = Dataset)) +
    geom_text(aes(label = Motif), position = position_quasirandom(), size = 3) +
    #geom_quasirandom(size = 0.5) +
    #geom_label_repel(aes(label = Motif), size = 2, label.size = 0.2, label.padding = 0.2) +
    facet_wrap(~Type, scales = "free", nrow = 1) +
    theme_bw()
  ggsave(paste0(data_path, "ngram_results/Most_frequent_motifs_", gsub(" ", "_", ith_type), ".png"), 
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


lapply(c("Pentagrams with gaps", "Tetragrams without gaps", 
         "Tetragrams with gaps", "Trigrams without gaps", "Trigrams with gaps"), function(ith_type) {
           lapply(unique(df_freq_decoded[["Encoding"]]), function(ith_enc) {
             x <- filter(df_freq_decoded, Type == ith_type, Encoding == ith_enc)
             p <- lapply(unique(x[["Dataset"]]), function(i) {
               filter(x, Dataset == i) %>% 
                 arrange(desc(Frequency)) %>% 
                 head(n = 15)
             }) %>% bind_rows() %>% 
               ggplot(aes(x = Frequency, y = reorder(Motif, Frequency))) +
               geom_col() +
               facet_wrap(~Dataset, scales = "free_y", nrow = 1) +
               theme_bw() +
               ggtitle(ith_type)
             ggsave(paste0(data_path, "ngram_results/Most_frequent_15best_", gsub(" ", "_", ith_type), "_", ith_enc, ".png"), 
                    p, width = 10, height = 5)
           })
         })



###--- Most frequent motifs, taxonomy ---###
df_freq_tax <- read_xlsx(paste0(data_path, "Motifs_results.xlsx"), sheet = "Taxonomy frequent")[,1:10] %>% 
  mutate(Motif = sapply(.[["Motif"]], function(i) gsub(".", " _ ", i, fixed = TRUE))) %>% 
  mutate(`Frequent in` = ifelse(`Frequent in` == "NA", "Unknown", `Frequent in`))

lapply(unique(df_freq_tax[["Type"]]), function(ith_type) {
  lapply(unique(df_freq_tax[["Dataset"]]), function(ith_set) {
    x <- filter(df_freq_tax, Dataset == ith_set, Type == ith_type)
    lapply(unique(x[["Frequent in"]]), function(ith_tax) {
      y <- filter(x, `Frequent in` == ith_tax)
      if(ith_tax %in% c("Viridiplantae", "Unknown") & ith_set == "cTP") {
        y <- setNames(y, c("Type", "Dataset", "Frequent in", "Motif", "Viridiplantae", "Unknown", "Freq3", "Freq4", "Frequency difference", "Cutoff"))
      } else if(ith_tax %in% c("Streptophyta", "Chlorophyta") & ith_set == "cTP") {
        y <- setNames(y, c("Type", "Dataset", "Frequent in", "Motif", "Streptophyta", "Chlorophyta", "Freq3", "Freq4", "Frequency difference", "Cutoff"))
      } else if(ith_set == "mTP") {
        y <- setNames(y, c("Type", "Dataset", "Frequent in", "Motif", "Viridiplantae", "Metazoa", "Fungi", "Unknown", "Frequency difference", "Cutoff"))
      } else if(ith_set == "SP") {
        y <- setNames(y, c("Type", "Dataset", "Frequent in", "Motif", "Eukaryota", "Bacteria", "Archaea", "Viruses", "Frequency difference", "Cutoff"))
      } 
      p <- ggplot(y, aes_string(x = ith_tax, y = paste0("reorder(Motif, ", ith_tax, ")"))) +
        geom_col() +
        facet_wrap(~Dataset, scales = "free_y", nrow = 1) +
        theme_bw() +
        ggtitle(paste0(ith_type, ", ", ith_set, ", ", ith_tax))
      ggsave(paste0(data_path, "ngram_results/Taxonomy_", gsub(" ", "_", ith_type), "_", ith_set, "_", ith_tax, ".png"), 
             p, width = 6, height = 2+nrow(y)*0.15, limitsize = FALSE)
    })
  })
})

lapply(unique(df_freq_tax[["Type"]]), function(ith_type) {
  lapply(unique(df_freq_tax[["Dataset"]]), function(ith_set) {
    x <- filter(df_freq_tax, Dataset == ith_set, Type == ith_type) %>% 
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
                                       Group == "Freq4" ~ "Unknown")) 
    } else if(ith_set == "SP") {
      y <- mutate(x, Group = case_when(Group == "Freq1" ~ "Eukaryota",
                                       Group == "Freq2" ~ "Bacteria",
                                       Group == "Freq3" ~ "Archaea",
                                       Group == "Freq4" ~ "Viruses")) 
    } 
    p <- ggplot(y, aes(y = Motif, x = Frequency, fill = Group)) +
      geom_col(position = position_dodge()) +
      facet_wrap(~`Frequent in`, scales = "free", nrow = 1) +
      theme_bw(base_size = 6) +
      ggtitle(paste0(ith_type, ", ", ith_set))
    ggsave(paste0(data_path, "ngram_results/Taxonomy_all_", gsub(" ", "_", ith_type), "_", ith_set, ".png"),
           p, width = 10, height = 2+nrow(y)*0.02, limitsize = FALSE)
  })
})




###--- Common motifs, taxonomy ---###
df_common_tax <- read_xlsx(paste0(data_path, "Motifs_results.xlsx"), sheet = "Taxonomy common") %>% 
  mutate(Motif = sapply(.[["Motif"]], function(i) gsub(".", " _ ", i, fixed = TRUE))) 

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
                                     Group == "Freq4" ~ "Unknown")) 
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
    ggtitle(ith_set)
  ggsave(paste0(data_path, "ngram_results/Taxonomy_common_", ith_set, ".png"),
         p, width = 4*length(unique(y[["Groups"]])), height = 4+nrow(y)*0.02, limitsize = FALSE)
})




###--- Differentiating motifs---###
df_diff <- read_xlsx(paste0(data_path, "Motifs_results.xlsx"), sheet = "Differences") %>% 
  mutate(Motif = sapply(.[["Motif"]], function(i) gsub(".", " _ ", i, fixed = TRUE))) 

lapply(unique(df_diff[["Type"]]), function(ith_type) {
  x <- filter(df_diff, Type == ith_type) %>% 
    select(-c("FreqDiff", "quipt_pval", "Cutoff")) %>% 
    pivot_longer(c(Freq1, Freq2), names_to = "Set", values_to = "Frequency")
  p <- ggplot(x, aes(y = Motif, x = Frequency, fill = Set)) +
    geom_col(position = position_dodge()) +
    facet_wrap(~Comparison, scales = "free", nrow = 1) +
    theme_bw() +
    ggtitle(ith_type)
  ggsave(paste0(data_path, "ngram_results/Differences_", ith_type, ".png"),
         p, width = 4*length(unique(x[["Comparison"]])), height = 4+nrow(x)*0.02, limitsize = FALSE)
})
