library(readr)
library(ggplot2)
library(scales)

df <- read_delim("impact_all_2.txt", delim = "\t", col_names = FALSE, show_col_types = FALSE)
colnames(df) <- c("Location", "Impact", "Consequence", "Type")

target_variants <- c("upstream_gene_variant", "intron_variant", "intergenic_variant", "downstream_gene_variant")

df$Impact_Grouped <- ifelse(df$Consequence %in% target_variants, 
                            df$Consequence, 
                            df$Impact)

df$Impact_Grouped <- factor(df$Impact_Grouped, 
                            levels = c("upstream_gene_variant", "downstream_gene_variant", 
                                       "intron_variant", "intergenic_variant", "MODIFIER",
                                       "HIGH", "MODERATE", "LOW"))

legend_mapping <- c(
  "mnps_biallelic_res"         = "MNP",
  "mnps_multiallelic_res"       = "MNP",
  "others_bigSV_res"           = "SV",
  "others_smallSV_res"         = "SCV",
  "others_snps_res"            = "SNV",
  "snps_biallelic_res"         = "SNV",
  "snps_multiallelic_res"       = "SNV",
  "indels_SV_private_res"      = "SV",
  "indels_small_SV_shared_res" = "SV",
  "indels_small_private_res"   = "InDel"
)

df$Type <- legend_mapping[df$Type]

p2 <- ggplot(df, aes(y = Impact_Grouped, fill = Type)) +
  geom_bar(position = "fill") +
  scale_x_continuous(labels = scales::percent) +
  labs(x = NULL, y = NULL, fill = NULL) + 
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(filename = "plot_percent_SV_effect.png", plot = p2, width = 6, height = 4, dpi = 300)
