#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=72:00:00
#SBATCH --mem=24GB

date

module purge
module use /apps/modules/all
module load BCFtools/1.17-GCC-11.2.0

dir="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v3/pggb"
out="ARVUO_pggb_all_chr_p95_s10k.vcf"
## concatenate the vcf files

#ls "$dir"/chr*/*.vcf.gz > vcf_list.txt ## make a list
bcftools concat -f vcf_list.txt -Oz -o "$out".gz ## merge

date
