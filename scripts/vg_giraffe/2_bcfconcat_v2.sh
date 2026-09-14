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

out_prefix="$1"
echo "Processing prefix: $out_prefix"

list_q10="${out_prefix}/tmp_q10_list.txt"
list_all="${out_prefix}/tmp_all_list.txt"

ls "$out_prefix"/*.q10.all.vcf.gz > "$list_q10" 2>/dev/null || true

> "$list_all"
for file in "$out_prefix"/*.all.vcf.gz; do
    if [ -f "$file" ] && [[ "$file" != *".q10."* ]]; then
        echo "$file" >> "$list_all"
    fi
done

if [ -s "$list_q10" ]; then
    echo "Running bcftools concat for Q10 files..."
    bcftools concat -f "$list_q10" -Oz -o "${out_prefix}/${out_prefix##*/}.q10.concat.vcf.gz"
else
    echo "No Q10 files found, skipping."
fi

if [ -s "$list_all" ]; then
    echo "Running bcftools concat for default ALL files..."
    bcftools concat -f "$list_all" -Oz -o "${out_prefix}/${out_prefix##*/}.default.concat.vcf.gz"
else
    echo "No default ALL files found, skipping."
fi

rm -f "$list_q10" "$list_all"

echo "Done!"
date
