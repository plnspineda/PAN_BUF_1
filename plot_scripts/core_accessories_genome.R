library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

setwd("/Users/polenpineda/Documents/buffalo_pangenome/dec_2025/core_accessories_genome/")
df <- "core_genomes_stats.txt"

df <- read_delim(df)

df_long_genome <- df %>%
  pivot_longer(cols = c(core, private_sum, shell_sum), names_to = "Type", values_to = "Value") %>%
  mutate(Type = factor(Type, levels = c("private_sum", "shell_sum", "core")))

df_swpc <- df %>%
  mutate(Type = "swpc_length") %>%
  select(chr, Type, Value = swpc_length)

df_rvuo <- df %>%
  mutate(Type = "rvuo_length") %>%
  select(chr, Type, Value = swpc_length)

# Plot
chr_levels <- c(as.character(1:23), "X")
plot_core_genomes <- ggplot() +
  geom_bar(data = df_long_genome, aes(x = factor(chr, levels = chr_levels), y = Value / 1e6, fill = factor(Type)),
           stat = "identity") +
  geom_bar(data = df_swpc, aes(x = factor(chr, levels = chr_levels), y = Value / 1e6, fill = Type),
           stat = "identity", position = position_nudge(x = 0.05), width = 0.3) +
  geom_bar(data = df_rvuo, aes(x = factor(chr, levels = chr_levels), y = Value / 1e6, fill = Type),
           stat = "identity", position = position_nudge(x = 0.30), width = 0.3) +
  labs(x = "Chromosome", y = "Size (Mb)") +
  scale_fill_manual(values = c("core" = "#427aa1", "shell_sum" = "#83c5be", "private_sum" = "#f4d35e", "swpc_length" = "#d1495b", "rvuo_length" = "darkblue"),
                    labels = c("core" = "Core", "shell_sum" = "Shell", "private_sum" = "Private", "swpc_length" = "PCC_UOA_SB_1v2", "rvuo_length" = "UOA_WB_1")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.title = element_blank()) +
  scale_y_continuous(labels = scales::comma, limits = c(0, 350))
plot_core_genomes
ggsave(plot_core_genomes, filename = "plot_core_genomes.png", width = 8, height = 4, dpi = 300)
