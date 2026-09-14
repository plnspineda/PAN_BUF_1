library(tidyverse)
library(ggplot2)
library(scales)

setwd("/Users/polenpineda/Documents/buffalo_pangenome/dec_2025/vep/vep_SVs/vep_files/txt_files/")
all_data_combined <- read_delim("all_data_vep_variant.clean.txt")

file_list <- c(
  "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.indels.SV.private.res.txt",
  "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.indels.small.SV.shared.res.txt",
  # "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.indels.small.private.res.txt",
  # "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.mnps.biallelic.res.txt",
  # "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.mnps.multiallelic.res.txt",
  "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.others.bigSV.res.txt"
  # "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.others.smallSV.res.txt",
  # "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.others.snps.res.txt",
  # "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.snps.biallelic.res.txt",
  # "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.snps.multiallelic.res.txt"
)

consequence_order <- c(
  "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
  "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
  "missense_variant", "inframe_insertion", "inframe_deletion",
  "protein_altering_variant", "splice_donor_5th_base_variant", "splice_region_variant",
  "splice_donor_region_variant", "splice_polypyrimidine_tract_variant",
  "start_retained_variant", "stop_retained_variant", "synonymous_variant",
  "coding_sequence_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
  "non_coding_transcript_exon_variant", "intron_variant", "non_coding_transcript_variant",
  "upstream_gene_variant", "downstream_gene_variant", "intergenic_variant"
)

target_variants <- c("upstream_gene_variant", "intron_variant", "intergenic_variant", "downstream_gene_variant")

legend_mapping <- c(
  "mnps.biallelic.res"         = "MNV",
  "mnps.multiallelic.res"       = "MNV",
  "others.bigSV.res"            = "SV",
  "others.smallSV.res"          = "SCV",
  "others.snps.res"             = "SNV",
  "snps.biallelic.res"          = "SNV",
  "snps.multiallelic.res"       = "SNV",
  "indels.SV.private.res"      = "SV",
  "indels.small.SV.shared.res" = "SV",
  "indels.small.private.res"   = "InDel"
)

### FUNCTIONS
process_vep_file <- function(file_path) {
  df <- read_tsv(
    file_path, comment = "#",
    col_names = c("Uploaded_variation", "Location", "Allele", "Gene", "Feature",
                  "Feature_type", "Consequence", "cDNA_position", "CDS_position",
                  "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra"),
    show_col_types = FALSE
  )

  unique_location_df <- df %>%
    mutate(First_Consequence = str_split_i(Consequence, ",", 1)) %>%
    mutate(severity_rank = match(First_Consequence, consequence_order)) %>%
    mutate(severity_rank = ifelse(is.na(severity_rank), Inf, severity_rank)) %>%
    group_by(Location) %>%
    slice_min(order_by = severity_rank, n = 1, with_ties = FALSE) %>%
    ungroup()

  processed_df <- unique_location_df %>%
    mutate(
      Impact = str_match(Extra, "IMPACT=([^;]+)")[, 2],
      Impact_Grouped = ifelse(First_Consequence %in% target_variants, First_Consequence, Impact),
      Impact_Grouped = factor(Impact_Grouped, levels = c(
        "upstream_gene_variant", "downstream_gene_variant",
        "intron_variant", "intergenic_variant", "MODIFIER",
        "HIGH", "MODERATE", "LOW"
      ))
    )

  return(processed_df)
}

SV_data_combined <- map_dfr(file_list, ~{
  file_key <- str_extract(.x, "(indels|mnps|others|snps)\\.[A-Za-z0-9_\\.]+\\.res")
  v_type <- legend_mapping[file_key]

  process_vep_file(.x) %>%
    mutate(
      File_Source = .x,
      Variant_Type = v_type
    )
})

summary_table <- all_data_combined %>% 
  group_by(Variant_Type, Impact_Grouped) %>% 
  summarise(count = n())

sum_impact <- summary_table %>% group_by(Impact_Grouped) %>%
  summarise(total = sum(count))

high_impact <- summary_table %>% filter(Impact_Grouped == "HIGH")

#write_tsv(all_data_combined, file = "all_data_vep_variant.clean.txt")

## plot the summary
summary_table$Impact_Grouped <- factor(summary_table$Impact_Grouped, levels = rev(c(
  "LOW", "MODERATE", "HIGH", "MODIFIER", 
  "intergenic_variant", "intron_variant", 
  "downstream_gene_variant", "upstream_gene_variant"
)))

my_colors <- c(
  "SV" = "#FF6B6B",
  "SNV"   = "#4D96FF",
  "SCV"   = "#6BCB40",
  "MNP"   = "#FFD93D",
  "InDel"    = "#B983FF")

p1 <- ggplot(summary_table, aes(y = Impact_Grouped, x = count, fill = Variant_Type)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) + 
  scale_x_continuous(labels = percent, expand = c(0.05, 0)) + 
  theme_classic() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal") +
  scale_fill_manual(values = my_colors)
p1

ggsave(p1, filename = "all_variant_VEP.png", width = 5, height = 3, dpi = 300)

### analyse the genes involve for high impact variants SVs

SV_high <- SV_data_combined %>% filter(Impact_Grouped == "HIGH")
SV_high_summary <- SV_high %>% group_by(Variant_Type, Impact_Grouped, File_Source) %>% 
  summarise(count = n())

### transcript ablations
### transcript ablations
gff <- read_delim("/Users/polenpineda/Documents/buffalo_pangenome/dec_2025/GO_enrich/gff_id.txt", col_names = c("Gene","Gene_Name"))
transcript_ablations <- SV_data_combined %>% filter(First_Consequence == "transcript_ablation")
transcript_ablations_2 <- merge(gff, transcript_ablations, by = "Gene", all.y = TRUE)
transcript_ablations_gene_counts <- transcript_ablations_2 %>%
  dplyr::mutate(Gene_Group = dplyr::if_else(stringr::str_detect(Gene_Name, "^LOC"), "LOC", Gene_Name)) %>%
  dplyr::count(Gene_Group, name = "Count")
