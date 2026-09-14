#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=24:00:00
#SBATCH --mem=16GB

date
module purge
module use /apps/modules/all
module load BCFtools/1.17-GCC-11.2.0
export TMPDIR=./

#bcftools view -f 'PASS' combined_snps.vcf.gz -O z -o combined_snps.pass.vcf.gz
#bcftools view -f 'PASS' combined_indels.vcf.gz -O z -o combined_indels.pass.vcf.gz
#bcftools view -f 'PASS' depth_combined_snps.vcf.gz -O z -o depth_combined_snps.pass.vcf.gz
bcftools view -f 'PASS' combined_snps.filter.vcf.gz -O z -o combined_snps.filter.pass.vcf.gz

date
