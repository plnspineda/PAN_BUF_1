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

echo "vcf_filtered"
bcftools concat -O z -o combined_snps.vcf.gz results/vcf_filtered/filtered_snps.{1..24}.vcf.gz results/vcf_filtered/filtered_snps.X.vcf.gz
bcftools concat -O z -o combined_indels.vcf.gz results/vcf_filtered/filtered_indels.{1..24}.vcf.gz results/vcf_filtered/filtered_indels.X.vcf.gz

echo "vcf_depth_filtered"
bcftools concat -O z -o depth_combined_snps.vcf.gz results/vcf_depth_filtered/depth_filtered_snps.{1..24}.vcf.gz results/vcf_depth_filtered/depth_filtered_snps.X.vcf.gz

date
