library(dplyr)
library(ggplot2)

data_path <- "/media/kasia/Data/Dropbox/Presequences/"

df_freq_tax2 <- readRDS("./data/df_freq_tax2.rds")


tax_data <- read.csv("./data/tax_plot_dat.csv")

class_org_numbers <- tax_data %>% 
  group_by(Dataset, class) %>% 
  summarise(count = length(unique(Organism))) %>% 
  setNames(c("Dataset", "tax", "count"))

phylum_org_numbers <- tax_data %>% 
  group_by(Dataset, phylum) %>% 
  summarise(count = length(unique(Organism))) %>% 
  setNames(c("Dataset", "tax", "count"))

kingdom_org_numbers <- tax_data %>% 
  group_by(Dataset, kingdom) %>% 
  summarise(count = length(unique(Organism))) %>% 
  setNames(c("Dataset", "tax", "count"))

superkingdom_org_numbers <- tax_data %>% 
  group_by(Dataset, superkingdom) %>% 
  summarise(count = length(unique(Organism))) %>% 
  setNames(c("Dataset", "tax", "count"))

all_org_numbers <- bind_rows(list(class_org_numbers,
                                  phylum_org_numbers,
                                  kingdom_org_numbers,
                                  superkingdom_org_numbers))

organism_motif_numbers <- lapply(unique(df_freq_tax2[["Type"]]), function(ith_type) {
  x <- filter(df_freq_tax2, Type == ith_type)
  lapply(unique(x[["Dataset"]]), function(ith_set) {
    y <- filter(x, Dataset == ith_set)
    lapply(unique(y[["Frequent in"]]), function(ith_tax) {
      z <- filter(y, `Frequent in` == ith_tax)
      data.frame(Type = ith_type,
                 Dataset = ith_set,
                 Group = ith_tax,
                 Motifs = length(unique(z[["Motif"]])),
                 Organisms = filter(all_org_numbers, Dataset == strsplit(ith_set, " ")[[1]][1], tax == ith_tax)[["count"]])
    }) %>% bind_rows()
  }) %>% bind_rows()
}) %>% bind_rows()


ggplot(organism_motif_numbers, aes(x = Organisms, y = Motifs, color = Dataset)) +
  geom_point() +
  theme_bw()

ggplot(organism_motif_numbers, aes(x = Organisms, y = Motifs, color = Dataset)) +
  geom_point() +
  theme_bw() + 
  facet_wrap(~Dataset) 

ggplot(organism_motif_numbers, aes(x = Organisms, y = Motifs, color = Group)) +
  geom_point() +
  theme_bw() + 
  facet_wrap(~Type)

ggplot(organism_motif_numbers, aes(x = Organisms, y = Motifs, color = Group)) +
  geom_point() +
  theme_bw() + 
  facet_wrap(~Dataset)

cor(organism_motif_numbers[["Motifs"]], organism_motif_numbers[["Organisms"]])
cor.test(organism_motif_numbers[["Motifs"]], organism_motif_numbers[["Organisms"]])

