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

dir="$1"
base=$(basename "$dir")
list_file="$dir/${base}_vcflist.txt"

echo "Compressing and indexing..."

find "$dir" -type f -name "*.vcf" -print0 | while IFS= read -r -d '' vcf; do
  echo "$vcf" >> "$list_file"
  echo "Compressing $vcf..."
  bgzip "$vcf"
  tabix -p vcf "${vcf}.gz"
done

date
