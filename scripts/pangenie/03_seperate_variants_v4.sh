#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=8:00:00
#SBATCH --mem=32GB

date
set -euo pipefail

module purge
module use /apps/modules/all
module load BCFtools/1.17-GCC-11.2.0

if [[ -z "${1:-}" ]]; then
    echo "Usage: $0 <vcf file>"
    exit 1
fi

vcf="$1"
base=$(basename "$vcf")
base=${base%%.vcf*}

echo "Input file: $vcf"
echo "Input prefix: $base"
echo "Output dir: $(pwd)"

if [[ ! -f "${vcf}.tbi" ]]; then
    tabix -p vcf -f "$vcf"
    echo "Index file: ${vcf}.tbi"
else
    echo "Index already exists, skipping: ${vcf}.tbi"
fi

### SNPs and MNPs
echo "Extract SNPs and MNPs"
bcftools view -v snps -Oz -o "${base}.snps.vcf.gz" "$vcf"
bcftools view -v mnps -Oz -o "${base}.mnps.vcf.gz" "$vcf"
tabix -p vcf "${base}.snps.vcf.gz"
tabix -p vcf "${base}.mnps.vcf.gz"

## SNPs with only 1bp ref and alt size
echo "Extract biallelic and multiallelic SNPs"
bcftools view -m2 -M2 -v snps -i 'strlen(REF)==1 && strlen(ALT)==1' -Oz -o "${base}.snps.biallelic.vcf.gz" "$vcf"
tabix -p vcf "${base}.snps.biallelic.vcf.gz"

bcftools view -v snps -i 'strlen(REF)==1 && strlen(ALT)==1' -m3 -Oz -o "${base}.snps.multiallelic.vcf.gz" "$vcf"
tabix -p vcf "${base}.snps.multiallelic.vcf.gz"

## MNPS with only 2bp ref and alt size
echo "Extract biallelic and multiallelic MNPs" ## ==2bp size
bcftools view -m2 -M2 -v mnps -i 'strlen(REF)==2 && strlen(ALT)==2' -Oz -o "${base}.mnps.biallelic.vcf.gz" "$vcf"
tabix -p vcf "${base}.mnps.biallelic.vcf.gz"

bcftools view -m3 -v mnps -i 'strlen(REF)==2 && strlen(ALT)==2' -Oz -o "${base}.mnps.multiallelic.vcf.gz" "$vcf"
tabix -p vcf "${base}.mnps.multiallelic.vcf.gz"

## InDels
echo "Extracting indels and other variants..."
bcftools view -v indels -Oz -o "${base}.indels.vcf.gz" "$vcf"
bcftools view -v other  -Oz -o "${base}.others.vcf.gz" "$vcf"
tabix -p vcf "${base}.indels.vcf.gz"
tabix -p vcf "${base}.others.vcf.gz"

echo "Extracting small indels..." ## get indels less than 50bp for both ref and alt
bcftools view -i 'strlen(REF)<50 && strlen(ALT)<50' -Oz -o "${base}.indels.small.vcf.gz" "${base}.indels.vcf.gz"
tabix -p vcf "${base}.indels.small.vcf.gz"

echo "Extracting SV indels..." ## get indels greater than 50bp for either ref or alt
bcftools view -i 'strlen(REF)>=50 || strlen(ALT)>=50' -Oz -o "${base}.indels.SV.vcf.gz" "${base}.indels.vcf.gz"
tabix -p vcf "${base}.indels.SV.vcf.gz"

##isec the indels to get non-overlapping variants (eg. variants that contains both small and SV indel)
bcftools isec -p isec_indels "${base}.indels.small.vcf.gz" "${base}.indels.SV.vcf.gz"
for f in isec_indels/*vcf; do echo -e "$(basename "$f")\t$(grep -vc '^#' "$f")"; done > isec_indels/count.txt
mv isec_indels/0000.vcf "${base}.indels.small.private.vcf"
mv isec_indels/0001.vcf "${base}.indels.SV.private.vcf"
mv isec_indels/0002.vcf "${base}.indels.small.SV.shared.vcf"

bgzip "${base}.indels.small.private.vcf"
bgzip "${base}.indels.SV.private.vcf"
bgzip "${base}.indels.small.SV.shared.vcf"

tabix -p vcf "${base}.indels.small.private.vcf.gz"
tabix -p vcf "${base}.indels.SV.private.vcf.gz"
tabix -p vcf "${base}.indels.small.SV.shared.vcf.gz"

bcftools concat "${base}.snps.biallelic.vcf.gz" "${base}.snps.multiallelic.vcf.gz" "${base}.mnps.biallelic.vcf.gz" "${base}.mnps.multiallelic.vcf.gz" "${base}.indels.small.private.vcf.gz" "${base}.indels.SV.private.vcf.gz" "${base}.indels.small.SV.shared.vcf.gz" -a -Oz -o "${base}.concat.vcf"
bgzip "${base}.concat.vcf"
tabix -p vcf -f "${base}.concat.vcf.gz"
bcftools isec -p isec_concat "$vcf" "${base}.concat.vcf.gz"
for f in isec_concat/*vcf; do echo -e "$(basename "$f")\t$(grep -vc '^#' "$f")"; done > isec_concat/count.txt
bgzip isec_concat/0000.vcf
tabix -p vcf -f isec_concat/0000.vcf.gz

bcftools view -v snps -Oz -o "${base}.others.snps.vcf.gz" isec_concat/0000.vcf.gz
tabix -p vcf "${base}.others.snps.vcf.gz"

bcftools isec -p isec_others_snps isec_concat/0000.vcf.gz "${base}.others.snps.vcf.gz"
mv isec_others_snps/0000.vcf "${base}.others.SVs.vcf"
bgzip "${base}.others.SVs.vcf"
tabix -p vcf -f "${base}.others.SVs.vcf.gz"

bcftools view -i 'strlen(REF)>=50 || strlen(ALT)>=50' -Oz -o "${base}.others.bigSV.vcf.gz" "${base}.others.SVs.vcf.gz"
tabix -p vcf -f "${base}.others.bigSV.vcf.gz"

bcftools isec -p isec_others_SVs "${base}.others.SVs.vcf.gz" "${base}.others.bigSV.vcf.gz"
mv isec_others_SVs/0000.vcf "${base}.others.smallSV.vcf"
bgzip "${base}.others.smallSV.vcf"
tabix -p vcf -f "${base}.others.smallSV.vcf.gz"

echo "Done!"
date
