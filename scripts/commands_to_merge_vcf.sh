conda activate pggb-env

# Pre-filter
for f in *.vcf; do
  bcftools view -i 'SVLEN!=0 && abs(SVLEN)>=50' "$f" > "${f%.vcf}.filt.vcf"
done

ls *.filt.vcf > vcf_list_all.txt

# Check SV types
for f in *.filt.vcf; do
	bcftools query -f '%SVTYPE\n' "$f"
done | sort | uniq -c

# 80723 DEL
#     1 DUP:INT
#   143 DUP:TANDEM
# 81409 INS

# Re-label DUP
for f in *.filt.vcf; do
  bcftools view "$f" | \
    awk 'BEGIN{OFS="\t"}
         /^#/ {print; next}
         {gsub(/SVTYPE=DUP:[A-Z]+/, "SVTYPE=DUP", $8); print}' \
  > "${f%.vcf}.norm.vcf"
done

ls *.filt.norm.vcf > vcf_list_all_norm.txt

# Check again the DUP has be relabelled properly
for f in *.filt.norm.vcf; do
        bcftools query -f '%SVTYPE\n' "$f"
done | sort | uniq -c

# 80723 DEL
#   144 DUP
# 81409 INS

conda deactivate
conda activate survivor

#survivor syntax
#SURVIVOR merge <list> <dist> <minSV> <type> <strand> <maxDist> <minSUP>

SURVIVOR merge vcf_list_all_norm.txt \
  500 \
  1 \
  1 \
  0 \
  50 \
  2 \
  merged_SVs.buffalo_pop.noINVBND.vcf

# Note: Chr is ARVUO#0#1 instead of 1. Renaming can be useful for some analyses.

conda deactivate
conda activate pggb-env

#sv length distribution
bcftools query -f '%INFO/SVTYPE\t%INFO/SVLEN\n' merged_SVs.buffalo_pop.noINVBND.vcf > sv_lengths.txt

# counts per SVTYPE
bcftools query -f '%INFO/SVTYPE\n' merged_SVs.buffalo_pop.noINVBND.vcf | sort | uniq -c
# 36279 DEL
#    90 DUP
# 44841 INS

# Standardize SVIM asm merged vcf to contain only INS and DEL
awk 'BEGIN{OFS="\t"}
/^##/ {print; next}
/^#CHROM/ {print $1,$2,$3,$4,$5,$6,$7,$8; next}
{
  if ($8 ~ /SVTYPE=DEL([;]|$)/ || $8 ~ /SVTYPE=INS([;]|$)/)
    print $1,$2,$3,$4,$5,$6,$7,$8
}' merged_SVs.buffalo_pop.noINVBND.vcf > merged_SVs.buffalo_pop.INSDEL.vcf

# SVIM merged
bcftools sort merged_SVs.buffalo_pop.INSDEL.vcf -o merged_SVs.buffalo_pop.INSDEL.sorted.vcf
bgzip merged_SVs.buffalo_pop.INSDEL.sorted.vcf
tabix -p vcf merged_SVs.buffalo_pop.INSDEL.sorted.vcf.gz

bcftools query -f '%CHROM\n' merged_SVs.buffalo_pop.INSDEL.sorted.vcf.gz | sort | uniq > svim_chroms.txt
################################################################################

#Bring in pggb + vcfbub file
# normalized multiallelic and keep biallelic
bcftools norm -m -any buffalo_pangenome_v3/SVIM/all_svim_vcf/matching_pggb/ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.vcf.gz | \
awk 'BEGIN{OFS="\t"}
/^#/ {print; next}
{
  ref = $4
  alt = $5
  ref_len = length(ref)
  alt_len = length(alt)

  # Classify SV type and calculate length change
  if (alt_len > ref_len) {
    svtype="INS"
    svlen = alt_len - ref_len
  } else if (ref_len > alt_len) {
    svtype="DEL"
    svlen = ref_len - alt_len
  } else {
    next # Skip SNPs or complex variants
  }

  # Size filter (>= 50 bp)
  if (svlen < 50) next

  # Annotate the INFO field cleanly
  if ($8 == "." || $8 == "") {
    $8 = "SVTYPE=" svtype ";SVLEN=" svlen
  } else {
    $8 = $8 ";SVTYPE=" svtype ";SVLEN=" svlen
  }

  print
}' > pggb_split_multiallelic_SVs.vcf
# Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:	27793942/1538313/0/0/0/0/0

# Extract the current header to a temporary file
bcftools view -h pggb_split_multiallelic_SVs.vcf > temp_header.txt

# Use sed to insert the new headers directly ABOVE the #CHROM line
sed -i '/^#CHROM/i ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">' temp_header.txt

# Apply the properly ordered header back onto your VCF data
bcftools reheader -h temp_header.txt pggb_split_multiallelic_SVs.vcf -o pggb_split_multiallelic_SVs.withheader.vcf

rm pggb_split_multiallelic_SVs.vcf

# PGGB / vcfbub
# sort and index
bcftools sort pggb_split_multiallelic_SVs.withheader.vcf -o pggb.sorted.vcf
bgzip pggb.sorted.vcf
tabix -p vcf pggb.sorted.vcf.gz
bcftools norm -m -both pggb.sorted.vcf.gz \
  -Oz -o pggb.norm.vcf.gz
tabix -p vcf pggb.norm.vcf.gz

# Using pggb as baseline for truvari
# Extract the list of chromosomes present in the PGGB VCF
bcftools query -f '%CHROM\n' pggb.norm.vcf.gz | sort | uniq > pggb_chroms.txt
# have a look at chr I have
# chr is matching svim, so below code is unnecessary
# zcat merged_SVs.buffalo_pop.INSDEL.sorted.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq | wc -l
#
# zcat merged_SVs.buffalo_pop.INSDEL.sorted.vcf.gz | grep -v "^#" | wc -l
#
#
# # Drop variants on unplaced, sex (X/Y), and MT chromosomes to get shared set
# bcftools view -e 'CHROM ~ "unplaced" || CHROM ~ "#X" || CHROM ~ "#Y" || CHROM ~ "#MT"' \
#   merged_SVs.buffalo_pop.INSDEL.sorted.vcf.gz \
#   -Oz -o merged_SVs.buffalo_pop.INSDEL.shared_chroms.vcf.gz
#
# tabix -p vcf merged_SVs.buffalo_pop.INSDEL.shared_chroms.vcf.gz
#
# zcat merged_SVs.buffalo_pop.INSDEL.shared_chroms.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq | wc -l
#
# zcat merged_SVs.buffalo_pop.INSDEL.shared_chroms.vcf.gz | grep -v "^#" | wc -l
#
#
# # Extract ONLY those dropped variants for separate backup file
# bcftools view -i 'CHROM ~ "unplaced" || CHROM ~ "#X" || CHROM ~ "#Y" || CHROM ~ "#MT"' \
#   merged_SVs.buffalo_pop.INSDEL.sorted.vcf.gz \
#   -Oz -o merged_SVs.buffalo_pop.INSDEL.unplaced_and_sex.vcf.gz
#
# tabix -p vcf merged_SVs.buffalo_pop.INSDEL.unplaced_and_sex.vcf.gz
#
# zcat merged_SVs.buffalo_pop.INSDEL.unplaced_and_sex.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq | wc -l
#
# zcat merged_SVs.buffalo_pop.INSDEL.unplaced_and_sex.vcf.gz | grep -v "^#" | wc -l

################################################################################
# truvari step
conda deactivate
conda activate truvari-env

# 1. Unfiltered PGGB as Baseline
truvari bench \
  -b pggb.norm.vcf.gz \
  -c merged_SVs.buffalo_pop.INSDEL.sorted.vcf.gz \
  -o truvari_pggb_base_unfiltered \
  --refdist 500 \
  --pctsize 0.5 \
  --pctseq 0.0 \
  --sizemin 50 \
  --passonly

# # 2. AF >= 0.05 PGGB as Baseline
# bcftools view -i 'AF>=0.05' pggb.norm.vcf.gz -Oz -o pggb.norm.AF0.05.vcf.gz
# tabix -p vcf pggb.norm.AF0.05.vcf.gz
#
# truvari bench \
#   -b pggb.norm.AF0.05.vcf.gz \
#   -c merged_SVs.buffalo_pop.INSDEL.shared_chroms.vcf.gz \
#   -o truvari_pggb_base_AF0.05 \
#   --refdist 500 \
#   --pctsize 0.5 \
#   --pctseq 0.0 \
#   --sizemin 50 \
#   --passonly
#
# # 3. AF >= 0.1 PGGB as Baseline
# bcftools view -i 'AF>=0.1' pggb.norm.vcf.gz -Oz -o pggb.norm.AF0.1.vcf.gz
# tabix -p vcf pggb.norm.AF0.1.vcf.gz
#
# truvari bench \
#   -b pggb.norm.AF0.1.vcf.gz \
#   -c merged_SVs.buffalo_pop.INSDEL.shared_chroms.vcf.gz \
#   -o truvari_pggb_base_AF0.1 \
#   --refdist 500 \
#   --pctsize 0.5 \
#   --pctseq 0.0 \
#   --sizemin 50 \
#   --passonly
#
# # 4. AF >= 0.2 PGGB as Baseline
# bcftools view -i 'AF>=0.2' pggb.norm.vcf.gz -Oz -o pggb.norm.AF0.2.vcf.gz
# tabix -p vcf pggb.norm.AF0.2.vcf.gz
#
# truvari bench \
#   -b pggb.norm.AF0.2.vcf.gz \
#   -c merged_SVs.buffalo_pop.INSDEL.shared_chroms.vcf.gz \
#   -o truvari_pggb_base_AF0.2 \
#   --refdist 500 \
#   --pctsize 0.5 \
#   --pctseq 0.0 \
#   --sizemin 50 \
#   --passonly
