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

echo "Extracting biallelic SNPs or InDels"
bcftools view -v snps,indels -m2 -M2 -Oz -o "${base}.biallelic.vcf.gz" "$vcf"
tabix -p vcf "${base}.biallelic.vcf.gz"

echo "Extracting multiallelic SNPs or InDels"
bcftools view -v snps,indels -m3 -Oz -o "${base}.multiallelic.vcf.gz" "$vcf"
tabix -p vcf "${base}.multiallelic.vcf.gz"

echo "Done!"
date
