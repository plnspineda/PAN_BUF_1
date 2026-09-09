library(dplyr)
library(ggplot2)
library(readr)
library(tidyverse)

setwd("/Users/polenpineda/Documents/buffalo_pangenome/dec_2025/")

df <- read_delim("alignment_summary_all_cleaned.csv")
fai <- read_delim("ARVUO.renamed.fa.fai", col_names = c("chromosome", "len", "coor", "X1", "X2"))
fai <- fai %>% select(chromosome, len)

df <- merge(df, fai, by = "chromosome")
df <- df %>% mutate(aligned_coverage = (Total_aligned*max_length)/len)
df <- df %>%
  mutate(chromosome = factor(chromosome,
                             levels = as.character(1:24)))

df_clean <- df %>%
  filter(!grepl("lib", sample))

df_per_chr <- df_clean %>%
  group_by(chromosome) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

plot_time_perchr <- df_clean %>%
  ggplot(aes(x= chromosome, y=Total_time_sec, group = sample, color = File)) +
  geom_line() +
  geom_point() +
  theme_classic()
plot_time_perchr

plot_time_perchr <- df_clean %>%
  ggplot(aes(x= chromosome, y=Speed_reads_per_sec, group = sample, color = File)) +
  geom_line() +
  geom_point() +
  theme_classic()
plot_time_perchr

plot_coverage_time <- df_clean %>%
  ggplot(aes(x=Total_time_sec, y=Total_alignments, color = File)) +
  geom_point() +
  theme_classic()
plot_coverage_time

df_clean_compute <- df_clean %>% mutate(aligned_percentage = Total_aligned/Total_alignments*100,
            properly_paired_percentage = Total_properly_paired/Total_paired*100)

df_summary <- df_clean %>% group_by(sample, File) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

df_summary_variants <- df_summary %>% select(File, sample, Insertions_bp, Deletions_bp, Substitutions_bp, Softclips_bp)
df_summary_variants <- df_summary %>% select(File, sample, Insertion_events, Deletion_events, Substitution_events, Softclip_events)

df_summary_variants_long <- df_summary_variants %>% pivot_longer(-c(File,sample), names_to = "metric", values_to = "value")

plot <- df_summary_variants_long %>%
  ggplot(aes(x = sample, y = value, fill = metric, color =File)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot
