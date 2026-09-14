#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=24:00:00
#SBATCH --mem=4GB

date
set -euo pipefail

module purge
module use /apps/modules/all
module load BCFtools/1.17-GCC-11.2.0

ls */*concat*vcf.gz > tmp.list

while read -r file; do
    tabix -p vcf "$file"
done < tmp.list

#find "$dir" -type f -name "*all.vcf.gz" > tmp.list
grep "default" tmp.list > default_vcf_list.txt
grep "q10" tmp.list > q10_vcf_list.txt

#echo "Compressing and indexing..."
echo "Run bcftools merge"
bcftools merge -l default_vcf_list.txt -Oz -o ARVUO_vggir_10x.default.all.vcf.gz
bcftools merge -l q10_vcf_list.txt -Oz -o ARVUO_vggir_10x.q10.all.vcf.gz
echo "Done!"

date
