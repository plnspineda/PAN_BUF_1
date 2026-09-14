#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time=72:00:00
#SBATCH --mem=64GB

date
module load BCFtools/1.17-GCC-11.2.0

for file in ARVUO_pangenie*.vcf.gz; do
    echo "Filtering $file..."
    bcftools view -O z -e 'AC==0 || AC==AN' "$file" -o "${file%.vcf.gz}.polymorphic.vcf.gz"
    bcftools index p_"${file%.vcf.gz}.polymorphic.vcf.gz"
done
date
