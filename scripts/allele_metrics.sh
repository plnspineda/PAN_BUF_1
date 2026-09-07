#!/usr/bin/env bash
set -euo pipefail

# =========================
# Output configuration
# =========================
OUT_TSV="pggb_all_chr_bub.allelicity_metrics.tsv"

# Print the header once at the start of the output file
echo -e "CHROM\tPOS\tALLELICITY\tALT_AF\tREF_AF\tAN\tVARIANT_CLASS" > "$OUT_TSV"

# =========================
# Process files in a loop
# =========================
# Loops over all vcf.gz files matching your pattern
for BUBBLE_VCF in ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.*.vcf.gz; do
    
    # Skip if no matching files are found
    [[ -e "$BUBBLE_VCF" ]] || continue

    # Extract the variant class from the filename
    # E.g., "ARVUO_..._r1M.indels.small.SV.shared.vcf.gz" -> "indels.small.SV.shared"
    vcf_basename=$(basename "$BUBBLE_VCF" .vcf.gz)
    variant_class=${vcf_basename#ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.}

    echo "[PROCESSING] File: $BUBBLE_VCF (Class: $variant_class)" >&2

    # Extract fields and stream directly into awk
    bcftools query \
      -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AN\n' \
      "$BUBBLE_VCF" \
    | awk -v v_class="$variant_class" '
    BEGIN {
        OFS = "\t";
    }
    {
        alts = $4;
        af_str = $5;
        AN = $6;

        # ------------------------
        # Split ALT alleles to count them
        # ------------------------
        n_alt = split(alts, alt, ",");

        # ------------------------
        # Allelicity determination
        # ------------------------
        if (n_alt == 1) {
            allelicity = "biallelic";
        } else if (n_alt > 1) {
            allelicity = "multiallelic";
        } else {
            allelicity = "monomorphic";
        }

        # ------------------------
        # Allele frequency calculations
        # ------------------------
        n_af = split(af_str, af_arr, ",");
        
        sum_af = 0;
        for (i = 1; i <= n_af; i++) {
            sum_af += af_arr[i];
        }
        
        ref_af = 1.0 - sum_af;
        if (ref_af < 0) ref_af = 0; 

        # ------------------------
        # Output results + variant class variable passed from bash
        # ------------------------
        print $1, $2, allelicity, af_str, ref_af, AN, v_class;
    }
    ' >> "$OUT_TSV"

done

echo "[DONE] All metrics combined and written to: $OUT_TSV"
