#!/usr/bin/env bash
set -euo pipefail

# =========================
# Input / output
# =========================
BUBBLE_VCF="ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.vcf.gz"
OUT_TSV="pggb_all_chr_bub.allelicity_variant_type_length_afreq.tsv"

# =========================
# Check input
# =========================
if [[ ! -f "$BUBBLE_VCF" ]]; then
    echo "[ERROR] Input VCF not found: $BUBBLE_VCF" >&2
    exit 1
fi

# =========================
# Extract minimal fields
# =========================
# We only need CHROM, POS, REF, ALT at the bubble level
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AN\n' \
  "$BUBBLE_VCF" \
| awk '
BEGIN {
    OFS = "\t";
    print "CHROM","POS","ALLELICITY","VARIANT_TYPE","VARIANT_LENGTH","ALT_AF","REF_AF","AN";
}
{
    ref = $3;
    alts = $4;
    ref_len = length(ref);

    # ------------------------
    # Split ALT alleles
    # ------------------------
    n_alt = split(alts, alt, ",");

    # ------------------------
    # Allelicity
    # ------------------------
    if (n_alt == 1)
        allelicity = "biallelic";
    else
        allelicity = "multiallelic";

    # ------------------------
    # SNP check
    # ------------------------
    is_snp = (ref_len == 1);
    for (i = 1; i <= n_alt; i++) {
        if (length(alt[i]) != 1)
            is_snp = 0;
    }

    # ------------------------
    # Length-based evaluation
    # ------------------------
    max_abs_diff = 0;
    has_ins = 0;
    has_del = 0;

    for (i = 1; i <= n_alt; i++) {
        alt_len = length(alt[i]);
        diff = alt_len - ref_len;

        if (diff >= 50) has_ins = 1;
        if (-diff >= 50) has_del = 1;

        abs_diff = (diff < 0) ? -diff : diff;
        if (abs_diff > max_abs_diff)
            max_abs_diff = abs_diff;
    }

    # ------------------------
    # Variant type (HPRC-style)
    # ------------------------
    if (is_snp) {
        variant_type = "SNP";
    }
    else if (max_abs_diff < 50) {
        variant_type = "small-indel";
    }
    else if (has_ins && has_del) {
        variant_type = "SV-OTHERS";
    }
    else if (has_ins) {
        variant_type = "SV-INS";
    }
    else {
        variant_type = "SV-DEL";
    }

    # ------------------------
    # Variant length per ALT
    # ------------------------
    len_str = "";
    for (i = 1; i <= n_alt; i++) {
        alt_len = length(alt[i]);
        diff = alt_len - ref_len;

        if (i == 1)
            len_str = diff;
        else
            len_str = len_str "," diff;
    }

    # ------------------------
    # Allele frequency
    # ------------------------
    af_str = $5;
    AN = $6;

    # ALT allele frequencies (preserve order)
    alt_af_str = af_str;

    # Compute REF allele frequency
    n_af = split(af_str, af_arr, ",");

    sum_af = 0;
    for (i = 1; i <= n_af; i++) {
        sum_af += af_arr[i];
    }

    ref_af = 1 - sum_af;

    print $1, $2, allelicity, variant_type, len_str, alt_af_str, ref_af, AN;
}
' > "$OUT_TSV"

echo "[DONE] Variant allelicity and type written to: $OUT_TSV"

