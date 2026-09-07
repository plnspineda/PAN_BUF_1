#!/bin/bash

date

module purge
module use /apps/modules/all
module load Singularity/3.10.5
module load BCFtools/1.17-GCC-11.2.0

vcf="$1"
vcf_out=$(basename "$vcf" .vcf)_vcfbub_L0_r1M.vcf
echo "Input file: $vcf"
echo "Output file: $vcf_out"

echo "Normalise vcf with vcfbub (remove overlapping variants)"

buffalo_pangenome/tools/vcfbub/vcfbub -l 0 -r 1000000 --input "$vcf" > "$vcf_out"
