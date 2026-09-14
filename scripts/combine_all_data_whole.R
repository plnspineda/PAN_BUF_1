library(tidyverse)

setwd("/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v3/vep_pggb/vep_files_chr8_p93/txt_files")

file_list <- c(
  "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.indels.SV.private.res.txt",
  "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.indels.small.SV.shared.res.txt",
  "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.indels.small.private.res.txt",
  "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.mnps.biallelic.res.txt",
  "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.mnps.multiallelic.res.txt",
  "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.others.bigSV.res.txt",
  "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.others.smallSV.res.txt",
  "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.others.snps.res.txt",
  "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.snps.biallelic.res.txt",
  "ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.snps.multiallelic.res.txt"
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
  "mnps.biallelic.res"         = "MNP",
  "mnps.multiallelic.res"       = "MNP",
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

all_data_combined <- map_dfr(file_list, ~{
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
  summarise(count = n(), .groups = "drop")

write_tsv(all_data_combined, file = "all_data_vep_variant.clean.txt")

