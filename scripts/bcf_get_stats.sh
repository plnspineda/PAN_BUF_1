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

for i in *vcf.gz; do
  echo "$i"
  base=$(basename "$i" .vcf.gz).stats
  bcftools stats -S - "$i" > "$base"
  echo "Done running bcftools stats -S - $i > $base"
done

date
