#!/bin/bash
#SBATCH -p a100cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=72:00:00
#SBATCH --mem=32GB

date
#conda activate odgi
module use /apps/modules/all
module load BEDTools/2.31.0-GCC-11.2.0
#set -euo pipefail

threads="8"

graph_og="$1"
fai="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v3/reference_fa/ARVUO_withX.fa.fai"
window="100000"
#step="90000"

echo "Input file: $graph_og"

## clean and prep inputs
chrom=$(odgi paths -i $graph_og -L | grep "ARVUO")
echo "$chrom"

chr="${chrom##*#}"
slide_out="slide_${chr}_100kwin.bed"
final_stats_out="odgi_${chr}_100kwin_stats.txt"

echo "Generating windows for $chrom..."
grep -w "$chrom" "$fai" | awk '{print $1 "\t" $2}' > "chrom_${chr}_size.txt"
bedtools makewindows -g "chrom_${chr}_size.txt" -w "$window" > "$slide_out"

rm "chrom_${chr}_size.txt"

> "$final_stats_out"

echo "Processing windows sequentially..."
window_count=0

while read -r chrom_name start end; do
    ((window_count++))

    tmp_bed=$(mktemp tmp_window_XXXXXX.bed)

    echo -e "${chrom_name}\t${start}\t${end}" > "$tmp_bed"
    tmp_og="tmp_extract_${chr}_${window_count}.og"

    echo "Processing row ${window_count}: ${chrom_name}:${start}-${end}"
    odgi extract -i "$graph_og" -b "$tmp_bed" -o "$tmp_og" -t "$threads" -E
    odgi_data=$(odgi stats -i "$tmp_og" | tail -n 1)
    echo -e "${chrom_name}:${start}-${end}\t${odgi_data}" >> "$final_stats_out"

    rm "$tmp_bed"
    rm "$tmp_og"

done < "$slide_out"

echo "Done! All stats are compiled in: $final_stats_out"
date
