rm(list = ls())

library(dplyr)
library(stringr)

setwd('C:/Users/tr432/Desktop/Roh/_Study_Design/25.05.02~_Gut-Lung_Axis/Exp1+1-S/MGX/04.metaphlan4.out')

df <- read.delim("merged_abundance_table_graphlan용.txt", check.names = FALSE)

sample_cols <- setdiff(colnames(df), c("Kingdom","Phylum","Class","Order","Family","Genus","Species","ASV"))

species_collapsed <- df %>%
  group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  summarise(across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")


write.table(species_collapsed, file = "species_collapsed.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
