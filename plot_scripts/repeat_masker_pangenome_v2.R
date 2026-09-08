library(dplyr)
library(readr)
library(stringr)
library(scales)
library(tidyr)
library(ggplot2)

setwd("/Users/polenpineda/Documents/buffalo_pangenome/dec_2025/repeat_masker/")
system("tail -n +3 SV_allindels.fa.out | grep -v '^$' | sed -e 's/  */\\t/g' | sed 's/^\\t//' > SV_allindels.clean.out")

repm <- "SV_allindels.clean.out"
repm <- read_tsv(repm, col_names = c("score", "divergence", "deletion", "insertion", "query_name", "query_begin", 
                                        "query_end", "query_left", "strand","repeat", "family", "repeat_begin", 
                                        "repeat_end","repeat_left","ID", "remarks"))
repm_df <- repm %>% separate(query_name, 
                                into = c("type","chr", "seq_start", "copy_ID","size", "copy_rep"), 
                                sep = "_", 
                                remove = FALSE) %>%
  mutate(chr = str_remove(chr, "ARVUO#0#"))

repm_df_filt <- repm_df %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  filter(perc_id > 0.6) %>%
  mutate(class = family) %>%
  separate(class, into = c("class","subclass"), sep = "/")

repm_df_subset <- repm_df_filt %>%
  dplyr::select(ID, type, chr, seq_start, size, `repeat`, family, query_align_len) %>%
  dplyr::mutate(
    seq_start = as.numeric(seq_start),
    size = as.numeric(size),
    seq_end = seq_start + size)

repm_df_subset$size_bin <- cut(
  repm_df_subset$size,
  breaks = c(50, 1000, 100000, 500000, Inf),
  labels = c("50-1,000", "1,000-100,000", "100,000-500,000", ">500,000"),
  include.lowest = TRUE)

## bubble plot of the indel sizes
plot_chr_indel <- repm_df_subset %>%
  ggplot(aes(x = seq_start/1e6, y = factor(chr, levels = as.character(1:23)), size = size_bin, color = type)) +
  geom_point(alpha = 0.3) +
  scale_size_manual(
    values = c("50-1,000" = 0.5, "1,000-100,000" = 1, "100,000-500,000" = 3, ">500,000" = 6),
    name = "Size Range") +
  scale_x_continuous(
    breaks = seq(0, 300, by = 50),
    labels = scales::comma_format(),
    limits = c(0, 300)) +
  labs(x = "Genomic Position (Mb)", y = "Chromosome", title = "Structural variations (InDel)") +
  theme_classic()
plot_chr_indel
ggsave(plot_chr_indel, file = "plot_chr_indel.tiff", width = 8, height = 4, dpi = 600)

## separate bubble plot del and ins

# Bubble plot for Deletions (DEL)
plot_chr_del <- repm_df_subset %>%
  filter(type == "del") %>%
  ggplot(aes(x = seq_start/1e6, y = factor(chr, levels = as.character(1:23)), size = size_bin)) +
  geom_point(alpha = 0.3, color = "red") +  # Color for deletions
  scale_size_manual(
    values = c("50-1,000" = 0.5, "1,000-100,000" = 1, "100,000-500,000" = 3, ">500,000" = 6),
    name = "Size Range"
  ) +
  scale_x_continuous(
    breaks = seq(0, 300, by = 50),
    labels = scales::comma_format(),
    limits = c(0, 300)
  ) +
  labs(x = "Genomic Position (Mb)", y = "Chromosome", title = "Structural Variations (Deletions)") +
  theme_classic()

# Bubble plot for Insertions (INS)
plot_chr_ins <- repm_df_subset %>%
  filter(type == "ins") %>%
  ggplot(aes(x = seq_start/1e6, y = factor(chr, levels = as.character(1:23)), size = size_bin)) +
  geom_point(alpha = 0.3, color = "blue") +  # Color for insertions
  scale_size_manual(
    values = c("50-1,000" = 0.5, "1,000-100,000" = 1, "100,000-500,000" = 3, ">500,000" = 6),
    name = "Size Range"
  ) +
  scale_x_continuous(
    breaks = seq(0, 300, by = 50),
    labels = scales::comma_format(),
    limits = c(0, 300)
  ) +
  labs(x = "Genomic Position (Mb)", y = "Chromosome", title = "Structural Variations (Insertions)") +
  theme_classic()

# Display the plots
plot_chr_del
plot_chr_ins

summary_type <- repm_df_filt %>%
  group_by(class) %>%
  summarise(bp = sum(query_align_len))

to_collate <- c("Low_complexity", "snRNA", "rRNA", "tRNA", "Unknown", "srpRNA", "RC", "scRNA")
summary_type_collated <- repm_df_filt %>%
  group_by(class) %>%
  summarise(bp = sum(query_align_len)) %>%
  mutate(class = if_else(class %in% to_collate, "Others", class)) %>%
  group_by(class) %>%
  summarise(bp = sum(bp)) %>%
  arrange(desc(bp))

repeat_summary <- repm_df_filt %>%
  group_by(`repeat`) %>%
  summarise(
    count = n(),
    total_size = sum(query_align_len, na.rm = TRUE),
    average_size = mean(query_align_len, na.rm = TRUE),
    max_size = max(query_align_len, na.rm = TRUE),
    min_size = min(query_align_len, na.rm = TRUE)) %>%
  arrange(desc(count))

# summary_type_edited <- data.frame(
#   class = c("LINE", "LTR", "SINE", "Simple_repeat", "Satellite", "DNA",
#             "Low_complexity", "Unknown", "RC", "RNA"),
#   bp = c(102499762, 19842672, 12224195, 4127481, 3310399, 1188579,
#          454924, 9345, 7696, 78286)
# )

plot_summary_type <- summary_type_collated %>%
  ggplot(aes(x = reorder(class, -bp), y = bp/1e6, fill = class)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(comma(bp/1e6, accuracy = 0.1))), 
            vjust = -0.5,
            size = 3.5) +
  theme_classic() +
  scale_y_continuous(labels = comma, 
                     expand = expansion(mult = c(0, 0.1))) +
  labs(y = "Size (Mbp)", x = "") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot_summary_type

ggsave(plot_summary_type, file = "plot_summary_type.png", width = 5, height = 3, dpi = 300)

#### filter for the highest peak

repm_df_filt_small <- repm_df_filt %>%
  filter(query_align_len > 120 & query_align_len < 150) %>%
  group_by(query_align_len, `repeat`) %>%
  summarise(count = n())

repm_df_filt_big <- repm_df_filt %>%
  filter(query_align_len > 1100 & query_align_len < 1400) %>%
  group_by(query_align_len, `repeat`) %>%
  summarise(count = n())

small_summ <- repm_df_filt_small %>%
  group_by(`repeat`) %>%
  summarise(total_count = sum(count), .groups = "drop")

big_summ <- repm_df_filt_big %>%
  group_by(`repeat`) %>%
  summarise(total_count = sum(count), .groups = "drop")

### plot the size of the SVs
## these SVs have higher count than the actual because I split multiallelic alleles (therefore, even if only one variant, if it's multiallelic, it can be two)
## the length is determined by subtracting the ref-alt, and then filtered to >50bp (since I have multiallelic data)

df_size <- read_delim("sv_sizes.txt",
                 delim = "\t",
                 col_names = c("chrom", "pos", "type", "len"),
                 show_col_types = FALSE)

plot_data <- df_size %>%
  filter(abs(len) >= 50) %>%
  mutate(
    size = abs(len),
    sv_type = if_else(len > 0, "Insertion", "Deletion")
  ) %>%
  group_by(size, sv_type) %>%
  summarise(count = n(), .groups = 'drop')

plot_size <- ggplot(plot_data, aes(x = size, y = count, color = sv_type)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  scale_x_log10(
    breaks = c(50, 100, 1000, 10000, 50000, 100000),
    labels = c("50 bp", "100 bp", "1 kbp", "10 kbp", "50 kbp", "100 kbp")) +
  scale_y_sqrt(
    breaks = c(0, 50, 100, 200, 400, 800, 1600)) +
  scale_color_manual(values = c("Deletion" = "#4A4A4A", "Insertion" = "#E69F00")) +
  labs(x = "Indel Size", y = "Count", color = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.85, 0.6),
        legend.text = element_text(size = 14),
        axis.title = element_text(face = "bold"))
plot_size
# 
# plot_data <- repm_df_filt %>%
#   mutate(size = as.numeric(size)) %>%
#   filter(size >= 50) %>%
#   group_by(size, type) %>%
#   summarise(count = n(), .groups = 'drop')
# 
# plot_size <- ggplot(plot_data, aes(x = size, y = count, color = type)) +
#   geom_line(linewidth = 0.6, alpha = 0.8) +
#   scale_x_log10(
#     breaks = c(50, 100, 1000, 10000, 50000, 100000),
#     labels = c("50 bp", "100 bp", "1 kbp", "10 kbp", "50 kbp", "100 kbp")
#   ) +
#   scale_color_brewer(palette = "Set1") +
#   labs(x = "Size", y = "Count", color = NULL) +
#   theme_classic(base_size = 14) +
#   theme(legend.position = c(0.85, 0.6),
#         legend.text = element_text(size = 14),
#         axis.title = element_text(face = "bold"))
# plot_size

ggsave(plot_size, filename="plot_size_indel_scaled.png", width = 8, height = 4, dpi = 300)
