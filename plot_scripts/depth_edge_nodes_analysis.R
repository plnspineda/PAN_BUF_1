library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(ggbreak)


setwd("/Users/polenpineda/Documents/buffalo_pangenome/dec_2025/complex_regions/")
df_compl <- read_delim("all_100k_stats.txt", col_names = c("chr","length","nodes","edges","paths","steps"))
df_depth <- read_delim("all_depth_ref_100kb_exact.bed")

bed_gaps <- read_delim("/Users/polenpineda/Documents/buffalo_pangenome/dec_2025/complex_regions/ARVUO_withX.gaps.coor",
                       col_names = c("chr_num","start","end"))

repmask <- "/Users/polenpineda/Documents/buffalo_pangenome/dec_2025/repeat_masker/UOA_WB_1_repeat.out.txt"
df_repmask <- read_delim(repmask, col_names = c("score", "divergence", "deletion", "insertion", "chr_num", "query_begin", 
                                                "query_end", "query_left", "strand","class", "family", "repeat_begin", 
                                                "repeat_end","repeat_left","ID", "remarks"))

fai <- read_delim("../ARVUO.fa.fai", col_names = c("chr","size"))
fai <- fai %>% select(chr, size)
chr_levels <- c(as.character(1:24), "X") 


df_depth <- df_depth %>%
  mutate(chr = `#path`)

df_compl <- df_compl %>%
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
df_compl <- df_compl %>%
  mutate(start = start +1, end = end +1)
df_compl <- merge(df_compl, fai, by = "chr")


df_compl_depth <- merge(df_compl, df_depth, by = c("chr", "start", "end"))
df_compl_depth <- df_compl_depth %>% distinct()

plot_compl_depth <- df_compl_depth %>%
  ggplot(aes(x = edge_node_ratio, y = mean.depth)) +
  geom_point(color = "steelblue", size = 2, alpha = 0.7) +
  scale_x_break(c(0.1, 1.15)) +
  scale_x_continuous(
    breaks = c(0.0, 1.2, 1.3, 1.4, 1.5),
    limits = c(0.0, 1.55)
  ) +
  theme_classic() +
  theme(
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    x = "Edge-to-node ratio",
    y = "Depth"
  )
plot_compl_depth

ggsave(plot_compl_depth, filename = "plot_compl_depth.png", dpi = 300, width = 4, height = 2.5)

### per chr
plot_compl_depth_chr <- df_compl_depth %>%
  ggplot(aes(x = edge_node_ratio, y = mean.depth)) +
  geom_point(color = "steelblue", size = 1, alpha = 0.7) +
  scale_x_continuous(
    breaks = c(0.0, 1.2, 1.3, 1.4, 1.5),
    limits = c(1.2, 1.55)
  ) +
  facet_wrap(~ chr_num) +
  theme_classic() +
  theme(
    axis.ticks = element_blank()
  ) +
  labs(
    x = "Edge-to-node ratio",
    y = "Depth"
  )

plot_compl_depth_chr

dgat1 <- df_compl_depth %>%
  filter(chr_num == 15, start > 81582793, end < 81793016)

### low mean depth but high edge-to-node ratio
#top10_depth_9_edge <- df_compl_depth %>%
#  arrange(abs(`mean.depth` - 9), desc(edge_node_ratio)) %>%
#  slice(1:10)

top_depth_edge <- df_compl_depth %>%
  filter(
    edge_node_ratio > 1.4,
    mean.depth > 9,
    mean.depth < 18
  ) %>% arrange(-edge_node_ratio) %>% slice(1:10)

#### repeats and complexity

df_depth_node <- df_compl_depth %>%
  select(chr_num, size, start, end, edge_node_ratio, depth = `mean.depth`) %>%
  mutate(chr_num = as.character(chr_num))

### clean repeatmasker result
df_repmask_clean <- df_repmask %>%
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  filter(perc_id > 0.6) %>%
  select(chr_num, start = query_begin, end = query_end,query_left,strand,class,family,query_align_len)


#### get intersect of the repeats and the depth range
library(fuzzyjoin)

bed_gaps_clean <- bed_gaps %>%
  mutate(chr_num = as.character(chr_num))

# Step 1: repmask (left) + depth (right), keep all repmask rows
df_merged <- genome_join(df_repmask_clean, df_depth_node,
                         by = c("chr_num", "start", "end"),
                         mode = "left")

df_merged <- df_merged %>%
  rename(
    chr_num    = chr_num.x,
    start      = start.x,
    end        = end.x,
    depth_start = start.y,
    depth_end   = end.y
  ) %>%
  select(-chr_num.y)   # duplicate key, not needed once .x is the reference

# Step 2: result (left) + bed_gaps_clean (right), keep all rows from step 1
df_merged <- genome_join(df_merged, bed_gaps_clean,
                         by = c("chr_num", "start", "end"),
                         mode = "left")

df_merged <- df_merged %>%
  rename(
    chr_num  = chr_num.x,
    start    = start.x,
    end      = end.x,
    gap_start = start.y,
    gap_end   = end.y
  ) %>%
  select(-chr_num.y)

head(df_merged)

summary_merged <- df_merged %>%
  group_by(chr_num, depth_start, depth_end, family) %>%
  summarise(length_repeats = sum(query_align_len), mean_depth = mean(depth),
            complexity = mean(edge_node_ratio), count = n())

summary_repeat_family <- df_merged %>%
  group_by(family) %>%
  summarise(
    length_repeats = sum(query_align_len, na.rm = TRUE),
    mean_depth = mean(depth, na.rm = TRUE),
    complexity = mean(edge_node_ratio, na.rm = TRUE),
    count = n()
  )

write_tsv(summary_repeat_family, "summary_repeat_family.tsv")


summary_length_repeats <- df_merged %>%
  group_by(chr_num, depth_start, depth_end) %>%
  summarise(length = sum(query_align_len),
            mean_depth = mean(depth),
            mean_complexity = mean(edge_node_ratio))

#mean(df_depth$mean.depth, na.rm = TRUE)
# 9.331415

## low edge to node ratio

low_edgenode <- df_merged %>%
  filter(chr_num == 17, depth_start == 72800001) 
# > sum(low_edgenode$query_align_len)
# [1] 54709
# > 54709/100000
# [1] 0.54709

no_edge <- df_merged %>%
  filter(edge_node_ratio < 1.0)
# > sum(no_edge$query_align_len)
# [1] 63749
# > 63749/100000*100
# [1] 63.749


chr8 <- df_repmask_clean %>% filter(chr_num == "8",
                                    start > 104845837, end < 105013602) #10kb flank
chr8_repeats <- chr8 %>% group_by(family) %>% summarise(length = sum(query_align_len))
# > 33479+1113+35121
# [1] 69713
# > 69713/147765
# [1] 0.4717829
# > 69713/147765*100
# [1] 47.17829
# > 101382-69713
# [1] 31669
# > 31669/147765*100
# [1] 21.432