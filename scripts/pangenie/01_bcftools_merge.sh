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

ls */*vcf > vcf_list.txt

echo "Decompressing and indexing..."
while read -r vcf; do
    bgzip -c "$vcf" > "${vcf}.gz"
    tabix -p vcf "${vcf}.gz"
done < vcf_list.txt

ls */*vcf.gz > vcfgz_list.txt
echo "Run bcftools merge"
bcftools merge -l vcfgz_list.txt -Oz -o ARVUO_pangenie_allsamples_withchr8_p93.vcf.gz

echo "Done!"
date
