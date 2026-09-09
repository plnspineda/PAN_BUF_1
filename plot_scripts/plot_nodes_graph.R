library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

setwd("/Users/polenpineda/Documents/buffalo_pangenome/dec_2025/complex_regions/")

df <- read_delim("all_100k_stats.txt", col_names = c("chr","length","nodes","edges","paths","steps"))
fai <- read_delim("../ARVUO.fa.fai", col_names = c("chr","size"))
fai <- fai %>% select(chr, size)

chr_levels <- c(as.character(1:24), "X") 
df_clean <- df %>%
  separate(chr, into = c("chr", "pos"), sep = ":") %>%
  separate(pos, into = c("start", "end"), sep = "-") %>%
  mutate(
    chr_num = sub("^ARVUO#0#", "", chr),
    chr_num = gsub("^(chr|Chr|chromosome)", "", chr_num),
    start = as.numeric(start),
    end = as.numeric(end),
    midpoint_mb = ((start + end) / 2) / 1e6,
    edge_node_ratio = edges / nodes
  ) %>%
  filter(chr_num %in% chr_levels) %>% 
  mutate(chr_num = factor(chr_num, levels = chr_levels))
df_clean <- merge(df_clean, fai, by = "chr")
## remove the chr8 outlier
df_clean <- df_clean %>% filter(edge_node_ratio > 1.0)

plot_violin_edge <- ggplot(df_clean, aes(x = chr_num, y = edge_node_ratio, fill = chr_num)) +
  geom_violin(alpha = 1, color = NA) +
  #geom_hline(yintercept = 1.364, color = "gray", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(width = 0.3, outlier.size = 0.1, fill = "white") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Chromosome",
    y = "Edge / Node Ratio"
  )
plot_violin_edge

##### top1%
chr_top1_thresholds <- df_clean %>%
  group_by(chr_num) %>%
  summarise(
    p99_threshold = quantile(edge_node_ratio, probs = 0.99, na.rm = TRUE),
    .groups = "drop"
  )

df_top1_per_chr <- df_clean %>%
  group_by(chr_num) %>%
  filter(edge_node_ratio >= quantile(edge_node_ratio, probs = 0.99, na.rm = TRUE)) %>%
  ungroup()


plot_violin_edge_99 <- ggplot(df_clean, aes(x = chr_num, y = edge_node_ratio, fill = chr_num)) +
  geom_violin(alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white") +
  geom_point(
    data = chr_top1_thresholds, 
    aes(x = chr_num, y = p99_threshold), 
    color = "maroon", shape = 17, size = 2
  ) +
  theme_classic() +
  scale_x_discrete(drop = FALSE) +
  theme(legend.position = "none") +
  labs(x = "Chromosome", y = "Edge / Node Ratio")
plot_violin_edge_99


#### 99% quartile that are at the terminal ends
df_top_terminal <- df_top1_per_chr %>%
  mutate(pos_end=size-end) %>%
  filter(end<15000000 | pos_end <15000000) %>%
  mutate(type = "terminal")

df_top_mid <- df_top1_per_chr %>%
  mutate(pos_end = size - end) %>%
  filter(start >= 15000000 & pos_end >= 15000000) %>%
  mutate(type = "mid")

df_top_label <- df_top1_per_chr %>%
  mutate(
    pos_end = size - end,
    type = if_else(start < 15e6 | pos_end < 15e6, "terminal", "mid")
  )

plot_node_to_edge <- ggplot(df_top_label, aes(x = midpoint_mb, y = chr_num, color = type)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("mid" = "#2b5c8f", "terminal" = "#e63946")) +
  theme_classic() +
  theme(
    # Adds horizontal grid lines
    panel.grid.major.y = element_line(color = "#e0e0e0", linewidth = 0.5),
    # Explicitly removes vertical lines (in case you switch base themes later)
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Genomic Position (Mb)",
    y = "Chromosome"
  )
plot_node_to_edge

ggsave(plot_node_to_edge, filename = "plot_node_to_edge_per_chr.png", dpi = 300, height = 3, width = 6)
ggsave(plot_violin_edge, filename = "plot_violin_edge.png", dpi = 300, height = 2.5, width = 4)
ggsave(plot_violin_edge_99, filename = "plot_violin_edge_99.png", dpi = 300, height = 3, width = 6)

write_tsv(df_clean, "df_clean.tsv")
getwd()


#####


