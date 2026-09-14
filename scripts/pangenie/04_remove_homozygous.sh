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

in="$1"
base=$(basename "$in" .vcf.gz)

bcftools view -i 'COUNT(GT="het")>0' "$in" -Oz -o "$base".nohomo.vcf.gz

echo "Done!"

date
